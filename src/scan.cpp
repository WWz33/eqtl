#include "eqtl/scan.hpp"
#include "eqtl/output.hpp"
#include "eqtl/stats_extra.hpp"
#include "eqtl/util.hpp"
#include "eqtl/plink_bed.hpp"
#include <fstream>
#include <random>
#include <atomic>
#include <future>
#include <sstream>
#include <mutex>
#include <thread>
#include <unordered_set>
#include <omp.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <cstring>
#include <utility>

namespace eqtl {

namespace {

struct ScopeOut {
  std::ofstream pairs;
  std::ofstream top;
  std::ofstream region;
  std::string tag;
  std::vector<char> pairs_buf;
  std::vector<char> top_buf;
  std::vector<char> region_buf;
};

bool needs_grm(Model m) {
  return m == Model::Lmm || m == Model::Glmm;
}

bool needs_counts(Model m) {
  return m == Model::Glm || m == Model::Glmm;
}

struct GeneReady {
  std::vector<int> keep;
  Eigen::VectorXd y;
  Eigen::MatrixXd X;
  Eigen::MatrixXd K;
  LmmBasis basis;
  bool has_basis = false;
};

bool build_gene_ready(const Eigen::VectorXd& y_full, const Eigen::MatrixXd& X_full,
                      const Eigen::MatrixXd* K_full, bool need_k, bool need_lmm_basis,
                      bool fast_sparse, GeneReady& out) {
  const int n_full = static_cast<int>(y_full.size());
  out.keep.clear();
  out.keep.reserve(static_cast<size_t>(n_full));
  for (int i = 0; i < n_full; ++i) {
    if (!std::isfinite(y_full(i))) continue;
    bool ok = true;
    for (int j = 0; j < X_full.cols(); ++j) {
      if (!std::isfinite(X_full(i, j))) {
        ok = false;
        break;
      }
    }
    if (ok) out.keep.push_back(i);
  }
  const int nk = static_cast<int>(out.keep.size());
  if (nk < 3) return false;

  out.y.resize(nk);
  out.X.resize(nk, X_full.cols());
  for (int r = 0; r < nk; ++r) {
    const int i = out.keep[static_cast<size_t>(r)];
    out.y(r) = y_full(i);
    out.X.row(r) = X_full.row(i);
  }
  if (out.y.array().abs().maxCoeff() < 1e-15) return false;

  if (need_k) {
    if (!K_full) return false;
    out.K.resize(nk, nk);
    for (int a = 0; a < nk; ++a)
      for (int b = 0; b < nk; ++b)
        out.K(a, b) = (*K_full)(out.keep[static_cast<size_t>(a)], out.keep[static_cast<size_t>(b)]);
    if (need_lmm_basis) {
      if (fast_sparse) {
        Eigen::MatrixXd K_use = out.K;
        sparsify_grm(K_use, 1e-4);
        out.basis = make_lmm_basis(K_use);
      } else {
        out.basis = make_lmm_basis(out.K);
      }
      out.has_basis = true;
    }
  }
  return true;
}

Eigen::VectorXd subset_dosage(const std::vector<double>& full, const std::vector<int>& keep) {
  Eigen::VectorXd g(static_cast<int>(keep.size()));
  for (size_t r = 0; r < keep.size(); ++r) {
    const int i = keep[r];
    if (i < 0 || static_cast<size_t>(i) >= full.size()) {
      g(static_cast<int>(r)) = std::numeric_limits<double>::quiet_NaN();
    } else {
      g(static_cast<int>(r)) = full[static_cast<size_t>(i)];
    }
  }
  return g;
}

AssocHit run_test(Model model, bool fast, const GeneReady& gr, const Eigen::VectorXd& g,
                  GenePrepLm* lm_cache, GenePrepLmm* lmm_cache, GenePrepGlm* glm_cache,
                  GenePrepGlmm* glmm_cache, bool have_cache) {
  switch (model) {
    case Model::Lm: {
      if (!have_cache) *lm_cache = prep_lm(gr.y, gr.X);
      return test_lm(*lm_cache, g);
    }
    case Model::Lmm: {
      if (!have_cache) {
        if (gr.has_basis) *lmm_cache = prep_lmm(gr.y, gr.X, gr.basis, fast);
        else *lmm_cache = prep_lmm(gr.y, gr.X, gr.K, fast);
      }
      return test_lmm(*lmm_cache, g);
    }
    case Model::Glm: {
      if (!have_cache) *glm_cache = prep_glm_nb(gr.y, gr.X, fast);
      return test_glm_nb(*glm_cache, g);
    }
    case Model::Glmm: {
      if (!have_cache) *glmm_cache = prep_glmm_pois(gr.y, gr.X, gr.K, fast);
      return test_glmm_pois(*glmm_cache, g);
    }
  }
  return {};
}

void prep_null(Model model, bool fast, const GeneReady& gr, GenePrepLm* lm_cache,
               GenePrepLmm* lmm_cache, GenePrepGlm* glm_cache, GenePrepGlmm* glmm_cache) {
  switch (model) {
    case Model::Lm:
      *lm_cache = prep_lm(gr.y, gr.X);
      break;
    case Model::Lmm:
      if (gr.has_basis) *lmm_cache = prep_lmm(gr.y, gr.X, gr.basis, fast);
      else *lmm_cache = prep_lmm(gr.y, gr.X, gr.K, fast);
      break;
    case Model::Glm:
      *glm_cache = prep_glm_nb(gr.y, gr.X, fast);
      break;
    case Model::Glmm:
      *glmm_cache = prep_glmm_pois(gr.y, gr.X, gr.K, fast);
      break;
  }
}

// Effect AF + variance on gene keep; monomorphic/non-finite → NaN.
// *maf_out = effect-allele frequency (mean dosage / 2); not folded to minor.
double subset_maf_or_nan(const Eigen::VectorXd& g, double* maf_out) {
  const int n = static_cast<int>(g.size());
  if (n <= 0) return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.0, sum2 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(g(i))) return std::numeric_limits<double>::quiet_NaN();
    sum += g(i);
    sum2 += g(i) * g(i);
  }
  double af = (sum / static_cast<double>(n)) / 2.0;
  if (af < 0) af = 0;
  if (af > 1) af = 1;
  if (maf_out) *maf_out = af;
  const double maf = af > 0.5 ? 1.0 - af : af;
  if (maf < 1e-12) return std::numeric_limits<double>::quiet_NaN();
  const double var = sum2 / n - (sum / n) * (sum / n);
  if (var < 1e-12) return std::numeric_limits<double>::quiet_NaN();
  return af;
}

// Test one SNP; fills meta only if p passes threshold (lazy string fill)
static void fill_snp_id(AssocHit& h, const SnpRec& snp) {
  if (!snp.id.empty()) h.snp = snp.id;
  else h.snp = snp.chrom + ":" + std::to_string(snp.pos) + ":" + snp.ref + ":" + snp.alt;
}

AssocHit test_one(Model model, bool fast, const GeneReady& gr, const SnpRec& snp,
                  const std::string& gene, const GeneLoc* loc, GenePrepLm* lm_c, GenePrepLmm* lmm_c,
                  GenePrepGlm* glm_c, GenePrepGlmm* glmm_c, bool have_cache,
                  Eigen::VectorXd& g_buf) {
  // reuse g_buf to avoid per-SNP alloc
  if (static_cast<int>(snp.dosage.size()) == static_cast<int>(gr.keep.size())) {
    if (g_buf.size() != static_cast<int>(snp.dosage.size()))
      g_buf.resize(static_cast<int>(snp.dosage.size()));
    std::memcpy(g_buf.data(), snp.dosage.data(), sizeof(double) * snp.dosage.size());
  } else {
    g_buf = subset_dosage(snp.dosage, gr.keep);
  }
  double maf_sub = snp.maf;
  if (!std::isfinite(subset_maf_or_nan(g_buf, &maf_sub))) {
    AssocHit h;
    h.p = std::numeric_limits<double>::quiet_NaN();
    return h;
  }
  AssocHit h = run_test(model, fast, gr, g_buf, lm_c, lmm_c, glm_c, glmm_c, have_cache);
  h.maf = maf_sub;
  h.n = static_cast<int>(gr.keep.size());
  // defer string fills — caller fills gene/snp/chrom/ref/alt only when needed
  h.gene = gene;
  h.chrom = snp.chrom;
  h.pos = snp.pos;
  h.ref = snp.ref;
  h.alt = snp.alt;
  if (loc && loc->ok) {
    h.has_tss_dist = true;
    h.tss_dist = static_cast<double>(snp.pos - loc->tss);
  }
  return h;
}

// Perm-only test: returns just p, no string alloc
double test_one_p(Model model, bool fast, const GeneReady& gr, const Eigen::VectorXd& g,
                  GenePrepLm* lm_c, GenePrepLmm* lmm_c, GenePrepGlm* glm_c, GenePrepGlmm* glmm_c) {
  double maf_sub;
  if (!std::isfinite(subset_maf_or_nan(g, &maf_sub))) return std::numeric_limits<double>::quiet_NaN();
  return run_test(model, fast, gr, g, lm_c, lmm_c, glm_c, glmm_c, true).p;
}

bool in_cis_window(const SnpRec& s, const GeneLoc& loc, int window) {
  if (!chrom_equal(s.chrom, loc.chrom)) return false;
  const int64_t cstart = std::max<int64_t>(1, loc.tss - window);
  const int64_t cend = loc.tss + window;
  return s.pos >= cstart && s.pos <= cend;
}


// Stream SNPs for one gene (list path removed; all callers stream).
template <typename StreamFn>
void scan_gene_snps(const Options& opt, Model model, const std::string& scope, const std::string& gene,
                    const GeneReady& gr, const GeneLoc* loc, double pthr, ScopeOut& out,
                    GeneSummary& summary,
                    StreamFn&& stream_snps) {
  GenePrepLm lm_c;
  GenePrepLmm lmm_c;
  GenePrepGlm glm_c;
  GenePrepGlmm glmm_c;
  prep_null(model, opt.fast, gr, &lm_c, &lmm_c, &glm_c, &glmm_c);

  AssocHit best;
  best.p = 2.0;
  AcatAcc acat_acc;
  int n_tested = 0;
  summary.n_sig = 0;
  Eigen::VectorXd g_buf;
  // perm: cache every tested dosage (exact min-p over full SNP set; O(M×n) when perm>0)
  std::vector<Eigen::VectorXd> cached_dosage;
  const bool do_perm = opt.perm > 0;

  auto apply_hit = [&](AssocHit& h, const SnpRec& snp) {
    if (!std::isfinite(h.p)) return;
    ++n_tested;
    acat_acc.add(h.p);
    if (h.p <= pthr) {
      fill_snp_id(h, snp);
      write_pair_line(out.pairs, h, model, scope);
      ++summary.n_sig;
    }
    if (h.p < best.p) {
      fill_snp_id(h, snp);
      best = h;
    }
  };

  stream_snps([&](const SnpRec& snp) {
    AssocHit h = test_one(model, opt.fast, gr, snp, gene, loc, &lm_c, &lmm_c, &glm_c, &glmm_c, true, g_buf);
    apply_hit(h, snp);
    if (do_perm && std::isfinite(h.p)) {
      cached_dosage.push_back(g_buf);
    }
  });

  summary.gene = gene;
  summary.n_tested = n_tested;
  if (loc && loc->ok) {
    summary.chrom = loc->chrom;
    summary.tss = loc->tss;
  }
  summary.acat_p = acat_acc.p();

  // Gene-level min-p perm (default --perm 0).
  // LM: residualize X then shuffle (Freedman–Lane).
  // GLM: residualize X, shuffle, re-add fit, round counts.
  // LMM: spectral null residual, shuffle, re-add null fit (K structure preserved better than OLS-on-y).
  // GLMM: residual y-mu_null, shuffle, re-add mu (approx).
  if (opt.perm > 0 && n_tested > 0) {
    const double T_obs = best.p;
    std::vector<double> T_perm(static_cast<size_t>(opt.perm));
    std::vector<double> perm_min_p(static_cast<size_t>(opt.perm));

    // LM/GLM: residualize fixed effects (Freedman–Lane style).
    // LMM: residualize in spectral space using null dinv weights (preserves K structure better than OLS).
    // GLMM: residualize on working scale y - mu_null (approx; still experimental).
    Eigen::VectorXd y_perm_base = gr.y;
    if (model == Model::Lm && lm_c.n > 0) {
      y_perm_base = lm_c.y_s;
    } else if (model == Model::Glm && gr.X.cols() > 0) {
      Eigen::LDLT<Eigen::MatrixXd> ldlt(gr.X.transpose() * gr.X);
      if (ldlt.info() == Eigen::Success) {
        const Eigen::VectorXd b = ldlt.solve(gr.X.transpose() * gr.y);
        y_perm_base = gr.y - gr.X * b;
      }
    } else if (model == Model::Lmm && lmm_c.n > 0) {
      // residualize y_til under null (X_til, dinv): y_til - X_til * beta0, back-transform by Q
      const Eigen::VectorXd& dinv = lmm_c.dinv;
      Eigen::MatrixXd XtDX = lmm_c.X_til.transpose() * dinv.asDiagonal() * lmm_c.X_til;
      Eigen::VectorXd XtDy = lmm_c.X_til.transpose() * (dinv.asDiagonal() * lmm_c.y_til);
      Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
      if (ldlt.info() == Eigen::Success) {
        const Eigen::VectorXd b0 = ldlt.solve(XtDy);
        const Eigen::VectorXd r_til = lmm_c.y_til - lmm_c.X_til * b0;
        y_perm_base = lmm_c.Q * r_til; // residual in original sample space
      }
    } else if (model == Model::Glmm && glmm_c.n > 0) {
      // Pearson-ish residual on mean scale: y - mu_null
      if (glmm_c.mu.size() == gr.y.size()) y_perm_base = gr.y - glmm_c.mu;
    }

    // ponytail: build grb once outside loop; only y changes per perm
    GeneReady grb;
    grb.keep = gr.keep;
    grb.X = gr.X;
    grb.K = gr.K;
    grb.basis = gr.basis;
    grb.has_basis = gr.has_basis;
    grb.y.resize(y_perm_base.size());

    std::atomic<int> perm_err{0};
#pragma omp parallel for schedule(dynamic) if (opt.threads > 1) firstprivate(grb)
    for (int b = 0; b < opt.perm; ++b) {
      if (perm_err.load()) continue;
      try {
      std::mt19937 rng_b(static_cast<unsigned>(opt.seed >= 0 ? opt.seed : 1) +
                         static_cast<unsigned>(b) * 9973u +
                         static_cast<unsigned>(std::hash<std::string>{}(gene)));

      std::vector<int> idx(static_cast<size_t>(y_perm_base.size()));
      std::iota(idx.begin(), idx.end(), 0);
      std::shuffle(idx.begin(), idx.end(), rng_b);

      for (int i = 0; i < y_perm_base.size(); ++i)
        grb.y(i) = y_perm_base(idx[static_cast<size_t>(i)]);
      if (model == Model::Glm && gr.X.cols() > 0) {
        // restore fixed-effect fit after residual shuffle; clamp to nonneg counts
        Eigen::LDLT<Eigen::MatrixXd> ldlt(gr.X.transpose() * gr.X);
        if (ldlt.info() == Eigen::Success) {
          const Eigen::VectorXd b = ldlt.solve(gr.X.transpose() * gr.y);
          const Eigen::VectorXd fit = gr.X * b;
          for (int i = 0; i < grb.y.size(); ++i) {
            double v = fit(i) + grb.y(i);
            if (v < 0) v = 0;
            grb.y(i) = std::round(v);
          }
        }
      } else if (model == Model::Glmm && glmm_c.mu.size() == grb.y.size()) {
        for (int i = 0; i < grb.y.size(); ++i) {
          double v = glmm_c.mu(i) + grb.y(i);
          if (v < 0) v = 0;
          grb.y(i) = std::round(v);
        }
      } else if (model == Model::Lmm && lmm_c.n > 0 && lmm_c.Q.size() > 0) {
        // residual was in sample space; add back null fixed fit in original y
        const Eigen::VectorXd& dinv = lmm_c.dinv;
        Eigen::MatrixXd XtDX = lmm_c.X_til.transpose() * dinv.asDiagonal() * lmm_c.X_til;
        Eigen::VectorXd XtDy = lmm_c.X_til.transpose() * (dinv.asDiagonal() * lmm_c.y_til);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
        if (ldlt.info() == Eigen::Success) {
          const Eigen::VectorXd b0 = ldlt.solve(XtDy);
          const Eigen::VectorXd fit = lmm_c.Q * (lmm_c.X_til * b0);
          grb.y += fit;
        }
      }

      GenePrepLm lm_b;
      GenePrepLmm lmm_b;
      GenePrepGlm glm_b;
      GenePrepGlmm glmm_b;
      prep_null(model, opt.fast, grb, &lm_b, &lmm_b, &glm_b, &glmm_b);

      double minp = 1.0;
      for (const auto& gd : cached_dosage) {
        const double p = test_one_p(model, opt.fast, grb, gd, &lm_b, &lmm_b, &glm_b, &glmm_b);
        if (std::isfinite(p) && p < minp) minp = p;
      }

      T_perm[static_cast<size_t>(b)] = -std::log10(std::max(minp, 1e-300));
      perm_min_p[static_cast<size_t>(b)] = minp;
      } catch (...) {
        perm_err.store(1);
      }
    }
    if (perm_err.load()) die("permutation failed (see earlier errors)");

    const double Tobs = -std::log10(std::max(T_obs, 1e-300));
    summary.p_emp = p_emp_count(Tobs, T_perm);
    if (!opt.disable_beta_approx) {
      beta_approx_p(perm_min_p, T_obs, summary.p_beta, summary.beta_shape1, summary.beta_shape2);
    } else {
      summary.p_beta = std::numeric_limits<double>::quiet_NaN();
    }
  } else {
    summary.p_emp = std::numeric_limits<double>::quiet_NaN();
    summary.p_beta = std::numeric_limits<double>::quiet_NaN();
  }

  if (best.p <= pthr && best.p <= 1.0) {
    write_pair_line(out.top, best, model, scope);
  }
}

// BH on a vector of GeneSummary: fills q_bh from acat_p
void bh_adjust(std::vector<GeneSummary>& gs) {
  const int G = static_cast<int>(gs.size());
  if (G == 0) return;
  std::vector<int> idx(static_cast<size_t>(G));
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return gs[a].acat_p < gs[b].acat_p; });
  double qmin = 1.0;
  for (int rank = G; rank >= 1; --rank) {
    const int i = idx[static_cast<size_t>(rank - 1)];
    double q = gs[i].acat_p * static_cast<double>(G) / static_cast<double>(rank);
    if (q > 1.0) q = 1.0;
    if (q < qmin) qmin = q;
    else q = qmin;
    gs[i].q_bh = q;
  }
}

}  // namespace


// SNP-outer LM for trans/gw (perm=0): one genotype pass. Same-keep → Y_s * g_s.
// Different keep / perm / non-LM → caller uses gene-outer.
struct GeneLmJob {
  std::string gene;
  GeneLoc loc;
  bool has_loc = false;
  GeneReady gr;
  GenePrepLm prep;
  GeneSummary summary;
  AssocHit best;
  AcatAcc acat_acc; // streaming — O(1) not O(M)
  // stage-2 perm (trans/gw, approximate): keep top-K dosages by nominal p
  std::vector<std::pair<double, Eigen::VectorXd>> top; // unsorted; trimmed to K
};

// max-heap by p (largest p at front) — retain smallest K p-values
static void topk_consider(std::vector<std::pair<double, Eigen::VectorXd>>& top, int K, double p,
                          const Eigen::VectorXd& g) {
  if (K <= 0 || !std::isfinite(p)) return;
  if (static_cast<int>(top.size()) < K) {
    top.emplace_back(p, g);
    if (static_cast<int>(top.size()) == K) {
      std::make_heap(top.begin(), top.end(),
                     [](const auto& a, const auto& b) { return a.first < b.first; });
    }
    return;
  }
  if (p >= top.front().first) return;
  std::pop_heap(top.begin(), top.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
  top.back() = {p, g};
  std::push_heap(top.begin(), top.end(),
                 [](const auto& a, const auto& b) { return a.first < b.first; });
}

template <typename Job>
static void fill_hit_meta(AssocHit& h, const Job& job, const SnpRec& snp) {
  h.gene = job.gene;
  h.chrom = snp.chrom;
  h.pos = snp.pos;
  h.ref = snp.ref;
  h.alt = snp.alt;
  fill_snp_id(h, snp);
  if (job.has_loc && job.loc.ok) {
    h.has_tss_dist = true;
    h.tss_dist = static_cast<double>(snp.pos - job.loc.tss);
  }
}

// per-job stats (thread-safe if each job owned by one thread); file write separate
template <typename Job>
static bool apply_snp_hit_stats(Job& job, AssocHit& h, const SnpRec& snp, double pthr) {
  if (!std::isfinite(h.p)) return false;
  ++job.summary.n_tested;
  job.acat_acc.add(h.p);
  bool write_pair = false;
  if (h.p <= pthr) {
    fill_hit_meta(h, job, snp);
    ++job.summary.n_sig;
    write_pair = true;
  }
  if (h.p < job.best.p) {
    fill_hit_meta(h, job, snp);
    job.best = h;
  }
  return write_pair;
}

template <typename Job>
static void apply_snp_hit(Job& job, AssocHit& h, const SnpRec& snp, double pthr, Model model,
                          ScopeOut& out) {
  if (apply_snp_hit_stats(job, h, snp, pthr)) {
    write_pair_line(out.pairs, h, model, out.tag);
  }
}

// Residualize g on prep.X (same as test_lm)
static bool residualize_g(const GenePrepLm& prep, const Eigen::VectorXd& g, Eigen::VectorXd& g_s,
                          double& gtg, Eigen::VectorXd& Xt_g_buf) {
  if (Xt_g_buf.size() != prep.p) Xt_g_buf.resize(prep.p);
  Xt_g_buf.noalias() = prep.X.transpose() * g;
  if (g_s.size() != g.size()) g_s.resize(g.size());
  g_s.noalias() = g - prep.X * (prep.XtX_inv * Xt_g_buf);
  gtg = g_s.squaredNorm();
  return gtg >= 1e-12;
}

static AssocHit hit_from_gty(const GenePrepLm& prep, double gtg, double gty, double maf_sub) {
  AssocHit h;
  h.n = prep.n;
  h.maf = maf_sub;
  if (gtg < 1e-12) {
    h.p = 1.0;
    return h;
  }
  h.beta = gty / gtg;
  const double df = prep.n - prep.p - 1;
  if (df <= 0) {
    h.p = 1.0;
    return h;
  }
  const double rss = prep.yty - h.beta * gty;
  const double s2 = std::max(rss / df, 0.0);
  h.se = std::sqrt(s2 / gtg);
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  h.r2 = (prep.yty > 0) ? (1.0 - rss / prep.yty) : 0.0;
  return h;
}

// Stage-2 gene perm on top-K SNPs only (trans/gw FastQTL-style; approximate, conservative).
// Document: p_emp uses min-p over top-K not full SNP set when --perm on trans/gw SNP-outer path.
template <typename Job>
void stage2_perm_topk(const Options& opt, Model model, Job& job,
                      const LmmBasis* ext_basis = nullptr) {
  if (opt.perm <= 0 || job.top.empty()) {
    job.summary.p_emp = std::numeric_limits<double>::quiet_NaN();
    job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  if (!(job.best.p < opt.perm_trans_thr) || job.summary.n_tested == 0) {
    job.summary.p_emp = std::numeric_limits<double>::quiet_NaN();
    job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  std::vector<Eigen::VectorXd> cached;
  cached.reserve(job.top.size());
  for (auto& e : job.top) cached.push_back(std::move(e.second));
  job.top.clear();

  GeneReady grb;
  grb.keep = job.gr.keep;
  grb.X = job.gr.X;
  grb.K = job.gr.K;
  if (ext_basis) {
    grb.basis = *ext_basis;
    grb.has_basis = true;
  } else {
    grb.basis = job.gr.basis;
    grb.has_basis = job.gr.has_basis;
  }
  grb.y = job.gr.y;

  GenePrepLm lm_c;
  GenePrepLmm lmm_c;
  GenePrepGlm glm_c;
  GenePrepGlmm glmm_c;
  prep_null(model, opt.fast, grb, &lm_c, &lmm_c, &glm_c, &glmm_c);

  Eigen::VectorXd y_perm_base = grb.y;
  if (model == Model::Lm && lm_c.n > 0) {
    y_perm_base = lm_c.y_s;
  } else if (model == Model::Lmm && lmm_c.n > 0) {
    const Eigen::VectorXd& dinv = lmm_c.dinv;
    Eigen::MatrixXd XtDX = lmm_c.X_til.transpose() * dinv.asDiagonal() * lmm_c.X_til;
    Eigen::VectorXd XtDy = lmm_c.X_til.transpose() * (dinv.asDiagonal() * lmm_c.y_til);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
    if (ldlt.info() == Eigen::Success) {
      const Eigen::VectorXd b0 = ldlt.solve(XtDy);
      const Eigen::VectorXd r_til = lmm_c.y_til - lmm_c.X_til * b0;
      y_perm_base = lmm_c.Q * r_til;
    }
  }

  const double T_obs = job.best.p;
  std::vector<double> T_perm(static_cast<size_t>(opt.perm));
  std::vector<double> perm_min_p(static_cast<size_t>(opt.perm));
  GeneReady grb2 = grb;
  grb2.y.resize(y_perm_base.size());
  std::atomic<int> perm_err{0};
#pragma omp parallel for schedule(dynamic) if (opt.threads > 1) firstprivate(grb2)
  for (int b = 0; b < opt.perm; ++b) {
    if (perm_err.load()) continue;
    try {
      std::mt19937 rng_b(static_cast<unsigned>(opt.seed >= 0 ? opt.seed : 1) +
                         static_cast<unsigned>(b) * 9973u +
                         static_cast<unsigned>(std::hash<std::string>{}(job.gene)));
      std::vector<int> idx(static_cast<size_t>(y_perm_base.size()));
      std::iota(idx.begin(), idx.end(), 0);
      std::shuffle(idx.begin(), idx.end(), rng_b);
      for (int i = 0; i < y_perm_base.size(); ++i)
        grb2.y(i) = y_perm_base(idx[static_cast<size_t>(i)]);
      if (model == Model::Lmm && lmm_c.n > 0 && lmm_c.Q.size() > 0) {
        const Eigen::VectorXd& dinv = lmm_c.dinv;
        Eigen::MatrixXd XtDX = lmm_c.X_til.transpose() * dinv.asDiagonal() * lmm_c.X_til;
        Eigen::VectorXd XtDy = lmm_c.X_til.transpose() * (dinv.asDiagonal() * lmm_c.y_til);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
        if (ldlt.info() == Eigen::Success) {
          const Eigen::VectorXd b0 = ldlt.solve(XtDy);
          grb2.y += lmm_c.Q * (lmm_c.X_til * b0);
        }
      }
      GenePrepLm lm_b;
      GenePrepLmm lmm_b;
      GenePrepGlm glm_b;
      GenePrepGlmm glmm_b;
      prep_null(model, opt.fast, grb2, &lm_b, &lmm_b, &glm_b, &glmm_b);
      double minp = 1.0;
      for (const auto& gd : cached) {
        const double p = test_one_p(model, opt.fast, grb2, gd, &lm_b, &lmm_b, &glm_b, &glmm_b);
        if (std::isfinite(p) && p < minp) minp = p;
      }
      T_perm[static_cast<size_t>(b)] = -std::log10(std::max(minp, 1e-300));
      perm_min_p[static_cast<size_t>(b)] = minp;
    } catch (...) {
      perm_err.store(1);
    }
  }
  if (perm_err.load()) die("stage-2 permutation failed");
  const double Tobs = -std::log10(std::max(T_obs, 1e-300));
  job.summary.p_emp = p_emp_count(Tobs, T_perm);
  if (!opt.disable_beta_approx) {
    beta_approx_p(perm_min_p, T_obs, job.summary.p_beta, job.summary.beta_shape1,
                  job.summary.beta_shape2);
  } else {
    job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
  }
}

template <typename G>
void scan_lm_snp_outer(const Options& opt, G& geno, const MissPolicy& mp, double maf,
                       const std::string& scope, ScopeOut& out, double pthr,
                       std::vector<GeneLmJob>& jobs) {
  if (jobs.empty()) return;

  // shared keep?
  bool same_keep = true;
  for (size_t i = 1; i < jobs.size(); ++i) {
    if (jobs[i].gr.keep != jobs[0].gr.keep) {
      same_keep = false;
      break;
    }
  }

  for (auto& j : jobs) {
    j.prep = prep_lm(j.gr.y, j.gr.X);
    j.best.p = 2.0;
    j.summary.gene = j.gene;
    j.summary.n_tested = 0;
    j.summary.n_sig = 0;
    if (j.has_loc && j.loc.ok) {
      j.summary.chrom = j.loc.chrom;
      j.summary.tss = j.loc.tss;
    }
  }

  if (same_keep) {
    // Y_s: G × n  (rows = residualized phenotypes)
    const int n = jobs[0].prep.n;
    const int Gz = static_cast<int>(jobs.size());
    Eigen::MatrixXd Ys(Gz, n);
    for (int gi = 0; gi < Gz; ++gi) Ys.row(gi) = jobs[static_cast<size_t>(gi)].prep.y_s.transpose();
    const auto& keep = jobs[0].gr.keep;
    const GenePrepLm& prep0 = jobs[0].prep;
    Eigen::VectorXd g_full(n), g_s(n), gty(Gz), Xt_g_buf(prep0.p);
    // ponytail: hoist once; per-SNP only fill/flag
    std::vector<AssocHit> write_hits(static_cast<size_t>(Gz));
    std::vector<char> write_flag(static_cast<size_t>(Gz), 0);
    const int topK = (opt.perm > 0) ? opt.perm_trans_top : 0;

    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      for (size_t r = 0; r < keep.size(); ++r) {
        const int i = keep[r];
        g_full(static_cast<int>(r)) =
            (i >= 0 && static_cast<size_t>(i) < snp.dosage.size())
                ? snp.dosage[static_cast<size_t>(i)]
                : std::numeric_limits<double>::quiet_NaN();
      }
      double maf_sub = snp.maf;
      if (!std::isfinite(subset_maf_or_nan(g_full, &maf_sub))) return true;

      double gtg = 0;
      if (!residualize_g(prep0, g_full, g_s, gtg, Xt_g_buf)) {
        for (int gi = 0; gi < Gz; ++gi) {
          auto& job = jobs[static_cast<size_t>(gi)];
          if (scope == "trans" && job.has_loc && in_cis_window(snp, job.loc, opt.window)) continue;
          AssocHit h;
          h.p = 1.0;
          h.n = prep0.n;
          h.maf = maf_sub;
          // p=1 monomorphic residual — skip top-K
          apply_snp_hit(job, h, snp, pthr, Model::Lm, out);
        }
        return true;
      }

      gty.noalias() = Ys * g_s;
      std::fill(write_flag.begin(), write_flag.end(), 0);
#pragma omp parallel for schedule(static) if (opt.threads > 1 && Gz > 32) num_threads(opt.threads)
      for (int gi = 0; gi < Gz; ++gi) {
        auto& job = jobs[static_cast<size_t>(gi)];
        if (scope == "trans" && job.has_loc && in_cis_window(snp, job.loc, opt.window)) continue;
        AssocHit h = hit_from_gty(job.prep, gtg, gty(gi), maf_sub);
        if (topK > 0 && std::isfinite(h.p) && h.p < 1.0) topk_consider(job.top, topK, h.p, g_full);
        if (apply_snp_hit_stats(job, h, snp, pthr)) {
          write_hits[static_cast<size_t>(gi)] = std::move(h);
          write_flag[static_cast<size_t>(gi)] = 1;
        }
      }
      for (int gi = 0; gi < Gz; ++gi) {
        if (write_flag[static_cast<size_t>(gi)])
          write_pair_line(out.pairs, write_hits[static_cast<size_t>(gi)], Model::Lm, out.tag);
      }
      return true;
    });
  } else {
    // different keep: still one I/O pass; reuse max-keep buffer
    size_t maxk = 0;
    for (const auto& j : jobs) maxk = std::max(maxk, j.gr.keep.size());
    Eigen::VectorXd g_buf(static_cast<int>(maxk));
    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      for (auto& job : jobs) {
        if (scope == "trans" && job.has_loc && in_cis_window(snp, job.loc, opt.window)) continue;
        const int nk = static_cast<int>(job.gr.keep.size());
        if (g_buf.size() < nk) g_buf.resize(nk);
        for (int r = 0; r < nk; ++r) {
          const int i = job.gr.keep[static_cast<size_t>(r)];
          g_buf(r) = (i >= 0 && static_cast<size_t>(i) < snp.dosage.size())
                         ? snp.dosage[static_cast<size_t>(i)]
                         : std::numeric_limits<double>::quiet_NaN();
        }
        Eigen::Map<Eigen::VectorXd> g(g_buf.data(), nk);
        double maf_sub = snp.maf;
        if (!std::isfinite(subset_maf_or_nan(g, &maf_sub))) continue;
        AssocHit h = test_lm(job.prep, g);
        h.maf = maf_sub;
        if (opt.perm > 0 && std::isfinite(h.p) && h.p < 1.0)
          topk_consider(job.top, opt.perm_trans_top, h.p,
                        Eigen::Map<const Eigen::VectorXd>(g_buf.data(), nk));
        apply_snp_hit(job, h, snp, pthr, Model::Lm, out);
      }
      return true;
    });
  }

  for (auto& job : jobs) {
    job.summary.acat_p = job.acat_acc.p();
    if (job.best.p <= pthr && job.best.p <= 1.0) {
      write_pair_line(out.top, job.best, Model::Lm, scope);
    }
    if (opt.perm > 0) stage2_perm_topk(opt, Model::Lm, job);
    else {
      job.summary.p_emp = std::numeric_limits<double>::quiet_NaN();
      job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
    }
  }
}


// SNP-outer LMM for trans/gw (perm=0). Same keep → share Q, g_til once per SNP.
// Different keep → still one genotype pass; per-gene test_lmm.
struct GeneLmmJob {
  std::string gene;
  GeneLoc loc;
  bool has_loc = false;
  GeneReady gr;
  GenePrepLmm prep;
  GeneSummary summary;
  AssocHit best;
  AcatAcc acat_acc;
  std::vector<std::pair<double, Eigen::VectorXd>> top;
};

// Workspace for LMM Wald — hoist out of per-gene / per-SNP loops
struct LmmTestWs {
  Eigen::MatrixXd Xg;
  Eigen::MatrixXd XtDX;
  Eigen::VectorXd XtDy, beta, e, cov_col;
};

// Wald on spectral space with precomputed g_til (same math as test_lmm)
static AssocHit test_lmm_gtil(const GenePrepLmm& prep, const Eigen::VectorXd& g_til, double maf_sub,
                              LmmTestWs& ws) {
  AssocHit h;
  h.n = prep.n;
  h.maf = maf_sub;
  const int p1 = prep.p + 1;
  if (ws.Xg.rows() != prep.n || ws.Xg.cols() != p1) ws.Xg.resize(prep.n, p1);
  ws.Xg.leftCols(prep.p) = prep.X_til;
  ws.Xg.col(prep.p) = g_til;
  const Eigen::VectorXd& dinv = prep.dinv;
  const double g_wss = g_til.dot(dinv.asDiagonal() * g_til);
  if (g_wss < 1e-12) {
    h.p = 1.0;
    return h;
  }
  if (ws.XtDX.rows() != p1 || ws.XtDX.cols() != p1) ws.XtDX.resize(p1, p1);
  if (ws.XtDy.size() != p1) ws.XtDy.resize(p1);
  ws.XtDX.noalias() = ws.Xg.transpose() * dinv.asDiagonal() * ws.Xg;
  ws.XtDy.noalias() = ws.Xg.transpose() * (dinv.asDiagonal() * prep.y_til);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(ws.XtDX);
  if (ldlt.info() != Eigen::Success) {
    h.p = 1.0;
    return h;
  }
  ws.beta = ldlt.solve(ws.XtDy);
  h.beta = ws.beta(prep.p);
  const int df = prep.n - prep.p - 1;
  double q = prep.y_til.dot(dinv.asDiagonal() * prep.y_til) - ws.XtDy.dot(ws.beta);
  if (q < 0) q = 0;
  const double sigma2 = (df > 0) ? (q / df) : 1.0;
  if (ws.e.size() != p1) ws.e.resize(p1);
  ws.e.setZero();
  ws.e(prep.p) = 1.0;
  ws.cov_col = ldlt.solve(ws.e);
  h.se = std::sqrt(std::max(sigma2 * ws.cov_col(prep.p), 0.0));
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  // partial R² (same as test_lmm)
  h.r2 = (prep.rss_null > 1e-15) ? std::max(0.0, 1.0 - q / prep.rss_null) : 0.0;
  return h;
}

// Free large GRM/basis after prep; keep y/X/keep for optional stage-2 perm
static void free_lmm_gene_raw(GeneLmmJob& j) {
  j.gr.K.resize(0, 0);
  j.gr.basis.Q.resize(0, 0);
  j.gr.basis.lambda.resize(0);
  j.gr.has_basis = false;
}

template <typename G>
void scan_lmm_snp_outer(const Options& opt, G& geno, const MissPolicy& mp, double maf,
                        const std::string& scope, ScopeOut& out, double pthr,
                        std::vector<GeneLmmJob>& jobs) {
  if (jobs.empty()) return;

  bool same_keep = true;
  for (size_t i = 1; i < jobs.size(); ++i) {
    if (jobs[i].gr.keep != jobs[0].gr.keep) {
      same_keep = false;
      break;
    }
  }

  LmmBasis shared_basis;
  bool have_shared_basis = false;
  if (same_keep && !jobs.empty()) {
    // one eigen-decomp for all genes with identical keep
    if (jobs[0].gr.has_basis) {
      shared_basis = jobs[0].gr.basis;
      have_shared_basis = true;
    } else if (jobs[0].gr.K.rows() > 0) {
      Eigen::MatrixXd K_use = jobs[0].gr.K;
      if (opt.fast) sparsify_grm(K_use, 1e-4);
      shared_basis = make_lmm_basis(K_use);
      have_shared_basis = true;
    }
  }

  for (auto& j : jobs) {
    if (have_shared_basis) j.prep = prep_lmm(j.gr.y, j.gr.X, shared_basis, opt.fast);
    else if (j.gr.has_basis) j.prep = prep_lmm(j.gr.y, j.gr.X, j.gr.basis, opt.fast);
    else j.prep = prep_lmm(j.gr.y, j.gr.X, j.gr.K, opt.fast);
    // stage-2 needs basis; free only when no perm
    if (opt.perm == 0) free_lmm_gene_raw(j);
    else {
      j.gr.K.resize(0, 0); // keep y/X; stage-2 uses ext_basis when same_keep
    }
    j.best.p = 2.0;
    j.summary.gene = j.gene;
    j.summary.n_tested = 0;
    j.summary.n_sig = 0;
    if (j.has_loc && j.loc.ok) {
      j.summary.chrom = j.loc.chrom;
      j.summary.tss = j.loc.tss;
    }
  }

  LmmTestWs ws;
  Eigen::VectorXd g_buf, g_til;

  if (same_keep) {
    for (size_t i = 1; i < jobs.size(); ++i) jobs[i].prep.Q.resize(0, 0);
    // keep shared_basis when perm>0 for stage-2 (~2MB); free Q copies on other jobs only
    if (opt.perm == 0) {
      shared_basis.Q.resize(0, 0);
      shared_basis.lambda.resize(0);
    }
    const Eigen::MatrixXd& Q = jobs[0].prep.Q;
    const auto& keep = jobs[0].gr.keep;
    const int n = static_cast<int>(keep.size());
    g_buf.resize(n);
    g_til.resize(n);
    const int Gz = static_cast<int>(jobs.size());
    std::vector<AssocHit> write_hits(static_cast<size_t>(Gz));
    std::vector<char> write_flag(static_cast<size_t>(Gz), 0);
    const int topK = (opt.perm > 0) ? opt.perm_trans_top : 0;

    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      for (size_t r = 0; r < keep.size(); ++r) {
        const int i = keep[r];
        g_buf(static_cast<int>(r)) =
            (i >= 0 && static_cast<size_t>(i) < snp.dosage.size())
                ? snp.dosage[static_cast<size_t>(i)]
                : std::numeric_limits<double>::quiet_NaN();
      }
      double maf_sub = snp.maf;
      if (!std::isfinite(subset_maf_or_nan(g_buf, &maf_sub))) return true;
      g_til.noalias() = Q.transpose() * g_buf;
      std::fill(write_flag.begin(), write_flag.end(), 0);
#pragma omp parallel if (opt.threads > 1 && jobs.size() > 32) num_threads(opt.threads)
      {
        LmmTestWs ws_t;
#pragma omp for schedule(static)
        for (int ji = 0; ji < Gz; ++ji) {
          auto& job = jobs[static_cast<size_t>(ji)];
          if (scope == "trans" && job.has_loc && in_cis_window(snp, job.loc, opt.window)) continue;
          AssocHit h = test_lmm_gtil(job.prep, g_til, maf_sub, ws_t);
          if (topK > 0 && std::isfinite(h.p) && h.p < 1.0) topk_consider(job.top, topK, h.p, g_buf);
          if (apply_snp_hit_stats(job, h, snp, pthr)) {
            write_hits[static_cast<size_t>(ji)] = std::move(h);
            write_flag[static_cast<size_t>(ji)] = 1;
          }
        }
      }
      for (int ji = 0; ji < Gz; ++ji) {
        if (write_flag[static_cast<size_t>(ji)])
          write_pair_line(out.pairs, write_hits[static_cast<size_t>(ji)], Model::Lmm, out.tag);
      }
      return true;
    });
    if (opt.perm == 0) jobs[0].prep.Q.resize(0, 0);
  } else {
    // one I/O; per-gene Q; reuse g_buf + g_til + LmmTestWs (no per-call Eigen alloc)
    size_t maxk = 0;
    for (const auto& j : jobs) maxk = std::max(maxk, j.gr.keep.size());
    g_buf.resize(static_cast<int>(maxk));
    g_til.resize(static_cast<int>(maxk));
    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      for (auto& job : jobs) {
        if (scope == "trans" && job.has_loc && in_cis_window(snp, job.loc, opt.window)) continue;
        const int nk = static_cast<int>(job.gr.keep.size());
        if (g_buf.size() < nk) g_buf.resize(nk);
        for (int r = 0; r < nk; ++r) {
          const int i = job.gr.keep[static_cast<size_t>(r)];
          g_buf(r) = (i >= 0 && static_cast<size_t>(i) < snp.dosage.size())
                         ? snp.dosage[static_cast<size_t>(i)]
                         : std::numeric_limits<double>::quiet_NaN();
        }
        double maf_sub = snp.maf;
        if (!std::isfinite(subset_maf_or_nan(
                Eigen::Map<const Eigen::VectorXd>(g_buf.data(), nk), &maf_sub)))
          continue;
        if (g_til.size() != nk) g_til.resize(nk);
        // g_til = Q' * g without temporary Map lifetime issues
        {
          Eigen::Map<const Eigen::VectorXd> g(g_buf.data(), nk);
          g_til.noalias() = job.prep.Q.transpose() * g;
        }
        AssocHit h = test_lmm_gtil(job.prep, g_til, maf_sub, ws);
        if (opt.perm > 0 && std::isfinite(h.p) && h.p < 1.0)
          topk_consider(job.top, opt.perm_trans_top, h.p,
                        Eigen::Map<const Eigen::VectorXd>(g_buf.data(), nk));
        apply_snp_hit(job, h, snp, pthr, Model::Lmm, out);
      }
      return true;
    });
  }

  for (auto& job : jobs) {
    job.summary.acat_p = job.acat_acc.p();
    if (job.best.p <= pthr && job.best.p <= 1.0) {
      write_pair_line(out.top, job.best, Model::Lmm, scope);
    }
    if (opt.perm > 0)
      stage2_perm_topk(opt, Model::Lmm, job, have_shared_basis ? &shared_basis : nullptr);
    else {
      job.summary.p_emp = std::numeric_limits<double>::quiet_NaN();
      job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
    }
  }
}


int run_make_grm(const Options& opt) {
  MissPolicy mp{opt.miss, opt.max_miss};
  const double maf = opt.maf;
  if (opt.use_bfile()) {
    PlinkBed bed;
    bed.open(opt.bfile);
    std::vector<std::string> samples = bed.samples();
    if (!opt.pheno.empty()) {
      PhenoData ph = load_pheno(opt.pheno);
      samples = intersect_order(ph.sample_ids, bed.samples());
      if (samples.empty()) die("no overlapping samples for GRM");
    }
    Grm g = compute_grm(bed, samples, mp, maf);
    write_grm_gcta(opt.out, g);
  } else {
    VcfSession vcf;
    vcf.open(opt.vcf);
    std::vector<std::string> samples = vcf.samples();
    if (!opt.pheno.empty()) {
      PhenoData ph = load_pheno(opt.pheno);
      samples = intersect_order(ph.sample_ids, vcf.samples());
      if (samples.empty()) die("no overlapping samples for GRM");
    }
    Grm g = compute_grm(vcf, samples, mp, maf);
    write_grm_gcta(opt.out, g);
  }
  return 0;
}


// cis: per-thread PlinkBed (independent FILE*) — genes sorted by chrom/TSS, chunked
static void run_cis_bfile_parallel(const Options& opt, PhenoData& ph, const CovData& cov,
                                   const std::unordered_map<std::string, GeneLoc>& annot,
                                   Eigen::MatrixXd* Kptr, bool need_k, bool need_lmm_basis,
                                   Model model, ScopeOut& so, double pthr,
                                   std::vector<GeneSummary>& summaries) {
  struct GeneWork {
    int gi = -1;
    std::string gene;
    GeneLoc loc;
  };
  std::vector<GeneWork> works;
  works.reserve(ph.gene_ids.size());
  for (size_t gi = 0; gi < ph.gene_ids.size(); ++gi) {
    const std::string& gene = ph.gene_ids[gi];
    auto it = annot.find(gene);
    if (it == annot.end()) continue;
    if (needs_counts(model) && !looks_like_counts(ph.Y.col(static_cast<int>(gi)))) {
      warn("skip non-count gene for " + model_str(model) + ": " + gene);
      continue;
    }
    works.push_back(GeneWork{static_cast<int>(gi), gene, it->second});
  }
  std::sort(works.begin(), works.end(), [](const GeneWork& a, const GeneWork& b) {
    if (a.loc.chrom != b.loc.chrom) return a.loc.chrom < b.loc.chrom;
    return a.loc.tss < b.loc.tss;
  });

  const int T = std::max(1, opt.threads);
  const MissPolicy mp{opt.miss, opt.max_miss};
  const double maf = opt.maf;
  std::vector<std::vector<GeneSummary>> part(static_cast<size_t>(T));
  std::vector<std::string> pairs_part(static_cast<size_t>(T));
  std::vector<std::string> top_part(static_cast<size_t>(T));
  std::atomic<int> err{0};

  LmmBasis grm_basis;
  bool have_grm_basis = false;
  if (need_lmm_basis && Kptr && Kptr->rows() == static_cast<int>(ph.sample_ids.size())) {
    Eigen::MatrixXd K_use = *Kptr;
    if (opt.fast) sparsify_grm(K_use, 1e-4);
    grm_basis = make_lmm_basis(K_use);
    have_grm_basis = true;
  }

#pragma omp parallel num_threads(T)
  {
    const int tid = omp_get_thread_num();
    const int nT = omp_get_num_threads();
    try {
      PlinkBed bed;
      bed.open(opt.bfile);
      bed.set_sample_order(ph.sample_ids);
      std::ostringstream pairs_ss, top_ss;
      // write headers only on merge; body only here
      ScopeOut local;
      // use stringstream via temporary files in memory — ScopeOut needs ostream
      // ponytail: thread-local ofstream to temp path
      const std::string t_pairs = opt.out + ".tmp." + std::to_string(tid) + ".pairs";
      const std::string t_top = opt.out + ".tmp." + std::to_string(tid) + ".top";
      const std::string t_region = opt.out + ".tmp." + std::to_string(tid) + ".region";
      local.tag = "cis";
      local.pairs.open(t_pairs);
      local.top.open(t_top);
      local.region.open(t_region);
      // no headers in temps

      const size_t nW = works.size();
      const size_t chunk = (nW + static_cast<size_t>(nT) - 1) / static_cast<size_t>(nT);
      const size_t lo = static_cast<size_t>(tid) * chunk;
      const size_t hi = std::min(nW, lo + chunk);
      for (size_t wi = lo; wi < hi; ++wi) {
        const auto& w = works[wi];
        Eigen::VectorXd y = ph.Y.col(w.gi);
        GeneReady gr;
        bool y_has_na = false;
        for (int i = 0; i < y.size(); ++i)
          if (!std::isfinite(y(i))) { y_has_na = true; break; }
        if (have_grm_basis && !y_has_na) {
          gr.keep.resize(static_cast<size_t>(y.size()));
          std::iota(gr.keep.begin(), gr.keep.end(), 0);
          gr.y = y;
          gr.X = cov.X;
          gr.basis = grm_basis;
          gr.has_basis = true;
        } else {
          if (!build_gene_ready(y, cov.X, Kptr, need_k, need_lmm_basis, opt.fast, gr)) continue;
        }
        GeneSummary summary;
        const GeneLoc* locp = &w.loc;
        const int64_t cstart = std::max<int64_t>(1, locp->tss - opt.window);
        const int64_t cend = locp->tss + opt.window;
        auto stream = [&](auto&& take) {
          bed.for_each_snp_region(locp->chrom, cstart, cend, mp, maf, [&](const SnpRec& s) {
            take(s);
            return true;
          });
        };
        scan_gene_snps(opt, model, "cis", w.gene, gr, locp, pthr, local, summary, stream);
        part[static_cast<size_t>(tid)].push_back(std::move(summary));
      }
      local.pairs.close();
      local.top.close();
      local.region.close();
      pairs_part[static_cast<size_t>(tid)] = t_pairs;
      top_part[static_cast<size_t>(tid)] = t_top;
      // drop empty region temps later
      (void)t_region;
    } catch (...) {
      err.store(1);
    }
  }
  if (err.load()) die("cis parallel scan failed");

  // merge pairs/top temps in thread order (= gene chrom/TSS order within chunks; chunks by sort)
  auto cat_file = [](std::ostream& dst, const std::string& path) {
    std::ifstream in(path);
    if (!in) return;
    dst << in.rdbuf();
    std::remove(path.c_str());
  };
  for (int tid = 0; tid < T; ++tid) {
    if (!pairs_part[static_cast<size_t>(tid)].empty())
      cat_file(so.pairs, pairs_part[static_cast<size_t>(tid)]);
    if (!top_part[static_cast<size_t>(tid)].empty())
      cat_file(so.top, top_part[static_cast<size_t>(tid)]);
    for (auto& s : part[static_cast<size_t>(tid)]) summaries.push_back(std::move(s));
    // cleanup region temp
    std::remove((opt.out + ".tmp." + std::to_string(tid) + ".region").c_str());
  }
}

template <typename G>
int run_eqtl_geno(const Options& opt, G& geno, PhenoData& ph,
                  const std::vector<std::string>& sample_order) {
  CovData cov = load_covar(opt.covar, sample_order);
  geno.set_sample_order(sample_order);

  std::unordered_map<std::string, GeneLoc> annot;
  const bool have_gff = !opt.gff.empty();
  if (have_gff) {
    annot = load_gff_tss(opt.gff, opt.gff_id_key);
    std::vector<std::string> gff_chroms;
    {
      std::unordered_set<std::string> seen;
      for (const auto& kv : annot) {
        if (seen.insert(kv.second.chrom).second)
          gff_chroms.push_back(kv.second.chrom);
      }
    }
    validate_chrom_names(geno.chromosomes(), gff_chroms);
  } else {
    info("no GFF: mode ignored; genome-wide all-pairs (gw)");
  }

  std::vector<std::string> scopes;
  if (!have_gff) {
    scopes = {"gw"};
  } else if (opt.mode == Mode::Cis) {
    scopes = {"cis"};
  } else if (opt.mode == Mode::Trans) {
    scopes = {"trans"};
  } else if (opt.mode == Mode::Gw) {
    scopes = {"gw"};
  } else {
    scopes = {"cis", "trans"};
  }

  bool any_grm = false;
  for (Model m : opt.models) {
    if (needs_grm(m)) any_grm = true;
  }
  Grm grm;
  Eigen::MatrixXd* Kptr = nullptr;
  Eigen::MatrixXd Kmat;
  MissPolicy mp{opt.miss, opt.max_miss};
  const double maf = opt.maf;
  if (any_grm) {
    if (!opt.grm.empty()) {
      grm = slice_grm(load_grm_gcta(opt.grm), sample_order);
    } else {
      info("computing GRM from genotypes (overlap samples)");
      grm = compute_grm(geno, sample_order, mp, maf);
    }
    Kmat = grm.K;
    Kptr = &Kmat;
  }

  // Never load_all into double matrix for trans/gw — stream instead.
  info(have_gff && scopes.size() == 1 && scopes[0] == "cis"
           ? "cis-only: stream region queries"
           : "trans/gw: stream SNPs (no load_all)");

  for (Model model : opt.models) {
    if (needs_grm(model) && !Kptr) die("internal: GRM required");
    const bool need_k = needs_grm(model);
    const bool need_lmm_basis = (model == Model::Lmm);

    for (const std::string& scope : scopes) {
      const std::string prefix = opt.out + "." + model_str(model) + "." + scope;
      ScopeOut so;
      so.tag = scope;
      so.pairs.open(prefix + ".pairs.tsv");
      so.top.open(prefix + ".top.tsv");
      so.region.open(prefix + ".region.tsv");
      if (!so.pairs || !so.top || !so.region) die("cannot open output for " + prefix);
      so.pairs_buf.assign(4 << 20, 0);
      so.top_buf.assign(1 << 16, 0);
      so.region_buf.assign(1 << 16, 0);
      so.pairs.rdbuf()->pubsetbuf(so.pairs_buf.data(), so.pairs_buf.size());
      so.top.rdbuf()->pubsetbuf(so.top_buf.data(), so.top_buf.size());
      so.region.rdbuf()->pubsetbuf(so.region_buf.data(), so.region_buf.size());
      write_pairs_header(so.pairs, model);
      write_top_header(so.top, model);
      write_region_header(so.region);

      const double pthr = (scope == "cis") ? opt.pval_cis : opt.pval_trans;
      std::vector<GeneSummary> summaries;

      // LM/LMM trans/gw: SNP-outer (I/O once). perm>0 → stage-2 top-K approx on genes with min-p < thr.
      const bool snp_outer =
          (scope == "trans" || scope == "gw") &&
          (model == Model::Lm || model == Model::Lmm);
      if (snp_outer && model == Model::Lm) {
        std::vector<GeneLmJob> jobs;
        jobs.reserve(ph.gene_ids.size());
        for (size_t gi = 0; gi < ph.gene_ids.size(); ++gi) {
          const std::string& gene = ph.gene_ids[gi];
          Eigen::VectorXd y = ph.Y.col(static_cast<int>(gi));
          GeneLoc loc_store;
          bool has_loc = false;
          if (have_gff) {
            auto it = annot.find(gene);
            if (it != annot.end()) {
              loc_store = it->second;
              has_loc = true;
            } else if (scope != "gw") {
              continue;
            }
          }
          if (scope == "trans" && !has_loc) continue;
          GeneLmJob job;
          job.gene = gene;
          job.loc = loc_store;
          job.has_loc = has_loc;
          if (!build_gene_ready(y, cov.X, Kptr, need_k, need_lmm_basis, opt.fast, job.gr)) continue;
          jobs.push_back(std::move(job));
        }
        info("trans/gw LM: SNP-outer (" + std::to_string(jobs.size()) + " genes)");
        scan_lm_snp_outer(opt, geno, mp, maf, scope, so, pthr, jobs);
        summaries.reserve(jobs.size());
        for (auto& j : jobs) summaries.push_back(std::move(j.summary));
      } else if (snp_outer && model == Model::Lmm) {
        std::vector<GeneLmmJob> jobs;
        jobs.reserve(ph.gene_ids.size());
        for (size_t gi = 0; gi < ph.gene_ids.size(); ++gi) {
          const std::string& gene = ph.gene_ids[gi];
          Eigen::VectorXd y = ph.Y.col(static_cast<int>(gi));
          GeneLoc loc_store;
          bool has_loc = false;
          if (have_gff) {
            auto it = annot.find(gene);
            if (it != annot.end()) {
              loc_store = it->second;
              has_loc = true;
            } else if (scope != "gw") {
              continue;
            }
          }
          if (scope == "trans" && !has_loc) continue;
          GeneLmmJob job;
          job.gene = gene;
          job.loc = loc_store;
          job.has_loc = has_loc;
          if (!build_gene_ready(y, cov.X, Kptr, need_k, need_lmm_basis, opt.fast, job.gr)) continue;
          jobs.push_back(std::move(job));
        }
        info("trans/gw LMM: SNP-outer (" + std::to_string(jobs.size()) + " genes)");
        scan_lmm_snp_outer(opt, geno, mp, maf, scope, so, pthr, jobs);
        summaries.reserve(jobs.size());
        for (auto& j : jobs) summaries.push_back(std::move(j.summary));
      } else if (scope == "cis" && opt.use_bfile() && opt.threads > 1 && have_gff) {
        info("cis: per-thread PlinkBed (" + std::to_string(opt.threads) + " threads)");
        run_cis_bfile_parallel(opt, ph, cov, annot, Kptr, need_k, need_lmm_basis, model, so, pthr,
                               summaries);
      } else {
        // shared GRM eigen once; only use when y has no NA (else build_gene_ready per gene)
        LmmBasis grm_basis;
        bool have_grm_basis = false;
        if (need_lmm_basis && Kptr && Kptr->rows() == static_cast<int>(ph.sample_ids.size())) {
          Eigen::MatrixXd K_use = *Kptr;
          if (opt.fast) sparsify_grm(K_use, 1e-4);
          grm_basis = make_lmm_basis(K_use);
          have_grm_basis = true;
          info("LMM: shared GRM eigen-decomp (n=" + std::to_string(Kptr->rows()) + ")");
        }

        // ponytail: no gene-level OMP — single FILE* bed + interleaved I/O; parallelize later per-thread bed
        for (size_t gi = 0; gi < ph.gene_ids.size(); ++gi) {
          const std::string& gene = ph.gene_ids[gi];
          Eigen::VectorXd y = ph.Y.col(static_cast<int>(gi));

          if (needs_counts(model) && !looks_like_counts(y)) {
            warn("skip non-count gene for " + model_str(model) + ": " + gene);
            continue;
          }

          const GeneLoc* locp = nullptr;
          GeneLoc loc_store;
          if (have_gff) {
            auto it = annot.find(gene);
            if (it != annot.end()) {
              loc_store = it->second;
              locp = &loc_store;
            } else if (scope != "gw") {
              continue;
            }
          }

          GeneReady gr;
          bool y_has_na = false;
          for (int i = 0; i < y.size(); ++i) {
            if (!std::isfinite(y(i))) { y_has_na = true; break; }
          }
          if (have_grm_basis && !y_has_na) {
            gr.keep.resize(static_cast<size_t>(y.size()));
            std::iota(gr.keep.begin(), gr.keep.end(), 0);
            gr.y = y;
            gr.X = cov.X;
            gr.basis = grm_basis;
            gr.has_basis = true;
          } else {
            if (!build_gene_ready(y, cov.X, Kptr, need_k, need_lmm_basis, opt.fast, gr)) continue;
          }

          summaries.emplace_back();
          GeneSummary& summary = summaries.back();

          if (scope == "cis") {
            if (!locp) { summaries.pop_back(); continue; }
            const int64_t cstart = std::max<int64_t>(1, locp->tss - opt.window);
            const int64_t cend = locp->tss + opt.window;
            auto stream = [&](auto&& take) {
              geno.for_each_snp_region(locp->chrom, cstart, cend, mp, maf, [&](const SnpRec& s) {
                take(s);
                return true;
              });
            };
            scan_gene_snps(opt, model, scope, gene, gr, locp, pthr, so, summary, stream);
          } else if (scope == "trans") {
            if (!locp) { summaries.pop_back(); continue; }
            auto stream = [&](auto&& take) {
              geno.for_each_snp(mp, maf, [&](const SnpRec& s) {
                if (in_cis_window(s, *locp, opt.window)) return true;
                take(s);
                return true;
              });
            };
            scan_gene_snps(opt, model, scope, gene, gr, locp, pthr, so, summary, stream);
          } else {
            auto stream = [&](auto&& take) {
              geno.for_each_snp(mp, maf, [&](const SnpRec& s) {
                take(s);
                return true;
              });
            };
            scan_gene_snps(opt, model, scope, gene, gr, locp, pthr, so, summary, stream);
          }
        }
      }
      bh_adjust(summaries);
      for (const auto& s : summaries) write_region_line(so.region, s);
      info("finished " + prefix);
    }
  }
  return 0;
}

int run_eqtl(const Options& opt) {
  omp_set_num_threads(opt.threads);

  PhenoData ph = load_pheno(opt.pheno);
  if (opt.use_bfile()) {
    PlinkBed bed;
    bed.open(opt.bfile);
    std::vector<std::string> sample_order = intersect_order(ph.sample_ids, bed.samples());
    if (sample_order.empty()) die("no overlapping samples between pheno and fam");
    info("overlap samples: " + std::to_string(sample_order.size()));
    slice_pheno(ph, sample_order);
    return run_eqtl_geno(opt, bed, ph, sample_order);
  }

  VcfSession vcf;
  vcf.open(opt.vcf);
  vcf.set_threads(opt.threads);
  std::vector<std::string> sample_order = intersect_order(ph.sample_ids, vcf.samples());
  if (sample_order.empty()) die("no overlapping samples between pheno and VCF");
  info("overlap samples: " + std::to_string(sample_order.size()));
  slice_pheno(ph, sample_order);
  return run_eqtl_geno(opt, vcf, ph, sample_order);
}

}  // namespace eqtl
