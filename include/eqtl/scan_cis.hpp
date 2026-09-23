/* eqtl — gene-outer cis scan with permutation (template on stream callback) */
#pragma once
#include "eqtl/scan_common.hpp"
#include <omp.h>
#include <random>
#include <atomic>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace eqtl {

// Gene-level cis scan: stream SNPs, test each, aggregate (ACAT), optional permutation.
// StreamFn signature: void(auto&& take) where take(const SnpRec&) is called per SNP.
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
  std::vector<Eigen::VectorXd> cached_dosage;
  // LMM keeps the genotype in the spectral domain instead. The permutation
  // scheme shuffles the phenotype's whitened residuals, never the genotype, so
  // g_til = Q^T g is the same on every draw — caching it removes one n²
  // projection per SNP per draw. maf_sub has to be stored alongside: it comes
  // from the sample-space dosage and cannot be recovered from g_til.
  std::vector<Eigen::VectorXd> cached_gtil;
  std::vector<double> cached_maf;
  const bool do_perm = opt.perm > 0;
  // Q is released on some trans paths before permutation; make sure the cis
  // projection below still has its full matrix.
  const bool cache_gtil =
      do_perm && model == Model::Lmm && lmm_c.n > 0 && lmm_c.Q.cols() == lmm_c.n;

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
      const size_t n_cached = cache_gtil ? cached_gtil.size() : cached_dosage.size();
      if (n_cached * static_cast<size_t>(std::max<int>(g_buf.size(), 1)) > 25'000'000)
        die("cis perm: too many SNPs in cis window for gene " + gene + " — reduce --window");
      if (cache_gtil) {
        cached_gtil.push_back(lmm_c.Q.transpose() * g_buf);
        cached_maf.push_back(h.maf);
      } else {
        cached_dosage.push_back(g_buf);
      }
    }
  });

  summary.gene = gene;
  summary.n_tested = n_tested;
  if (loc && loc->ok) {
    summary.chrom = loc->chrom;
    summary.tss = loc->tss;
  }
  summary.acat_p = acat_acc.p();

  // Gene-level min-p permutation.
  if (opt.perm > 0 && n_tested > 0) {
    const double T_obs = best.p;
    std::vector<double> T_perm(static_cast<size_t>(opt.perm));
    std::vector<double> perm_min_p(static_cast<size_t>(opt.perm));

    Eigen::VectorXd y_perm_base = gr.y;
    Eigen::VectorXd Xb0_til;  // LMM: null-fitted spectral mean X_til·b0
    if (model == Model::Lm && lm_c.n > 0) {
      y_perm_base = lm_c.y_s;
    } else if (model == Model::Lmm && lmm_c.n > 0) {
      // whitened spectral residuals: Var(w_i) = σ² under null → exchangeable.
      // b0 = A00^{-1} X_til^T D y_til does not depend on the permutation and is
      // already cached at prep as chi0. Both this block and the one inside the
      // loop below used to re-form and re-factor A00 — here once per gene, and
      // there again on every permutation draw.
      const Eigen::VectorXd& dinv = lmm_c.dinv;
      if (lmm_c.has_a00) {
        Xb0_til = lmm_c.X_til * lmm_c.chi0;
      } else {
        Eigen::MatrixXd XtDX = lmm_c.X_til.transpose() * dinv.asDiagonal() * lmm_c.X_til;
        Eigen::VectorXd XtDy = lmm_c.X_til.transpose() * (dinv.asDiagonal() * lmm_c.y_til);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
        if (ldlt.info() == Eigen::Success) Xb0_til = lmm_c.X_til * ldlt.solve(XtDy);
      }
      if (Xb0_til.size() == dinv.size()) {
        const Eigen::VectorXd r_til = lmm_c.y_til - Xb0_til;
        y_perm_base = r_til.cwiseProduct(dinv.cwiseSqrt());
      }
    } else if (model == Model::Glmm && glmm_c.n > 0) {
      if (glmm_c.mu.size() == gr.y.size()) y_perm_base = gr.y - glmm_c.mu;
    }
    // Glm: parametric NB sampling below; y_perm_base unused.

    GeneReady grb;
    grb.keep = gr.keep;
    grb.X = gr.X;
    grb.has_basis = gr.has_basis;
    // Follow the source's shared refs; gr.K / gr.basis are EMPTY on the shared-GRM
    // fast path, so &gr.K would dangle into a zero-size matrix for GLMM PQL.
    grb.K_ref = gr.K_ref ? gr.K_ref : &gr.K;
    grb.basis_ref = gr.basis_ref ? gr.basis_ref : &gr.basis;
    grb.y.resize(y_perm_base.size());

    std::atomic<int> perm_err{0};
    LmmTestWs ws_b;  // per-thread copy (firstprivate); see the trans SNP-outer loop
#pragma omp parallel for schedule(dynamic) if (opt.threads > 1 && !omp_in_parallel()) firstprivate(grb, ws_b)
    for (int b = 0; b < opt.perm; ++b) {
      if (perm_err.load()) continue;
      try {
        std::mt19937 rng_b(static_cast<unsigned>(opt.seed >= 0 ? opt.seed : 1) +
                           static_cast<unsigned>(b) * 9973u + fnv1a(gene));

        std::vector<int> idx(static_cast<size_t>(y_perm_base.size()));
        std::iota(idx.begin(), idx.end(), 0);
        std::shuffle(idx.begin(), idx.end(), rng_b);

        for (int i = 0; i < y_perm_base.size(); ++i)
          grb.y(i) = y_perm_base(idx[static_cast<size_t>(i)]);
        if (model == Model::Glm && glm_c.n > 0 && glm_c.mu.size() == grb.y.size()) {
          // parametric permutation: y* ~ NB(mu_null, phi_null)
          for (int i = 0; i < grb.y.size(); ++i) {
            const double mu_i = glm_c.mu(i);
            if (mu_i <= 1e-12) { grb.y(i) = 0; continue; }
            if (glm_c.phi > 1e-9) {
              std::gamma_distribution<double> gam(1.0 / glm_c.phi, glm_c.phi * mu_i);
              std::poisson_distribution<long> pois(gam(rng_b));
              grb.y(i) = static_cast<double>(pois(rng_b));
            } else {
              std::poisson_distribution<long> pois(mu_i);
              grb.y(i) = static_cast<double>(pois(rng_b));
            }
          }
        } else if (model == Model::Glmm && glmm_c.mu.size() == grb.y.size()) {
          for (int i = 0; i < grb.y.size(); ++i) {
            double v = glmm_c.mu(i) + grb.y(i);
            if (v < 0) v = 0;
            grb.y(i) = std::round(v);
          }
        } else if (model == Model::Lmm && lmm_c.n > 0 && lmm_c.Q.size() > 0) {
          // Xb0_til is the null-fitted spectral mean hoisted above. The size
          // guard covers it; the "no value" state cannot reach here, since it
          // only arises for a non-ok prep, and that yields no finite p, so the
          // permutation block (n_tested > 0) is never entered.
          if (Xb0_til.size() == lmm_c.dinv.size()) {
            const Eigen::VectorXd& dinv = lmm_c.dinv;
            // grb.y holds shuffled whitened residuals w_perm; un-whiten → spectral residual
            Eigen::VectorXd r_til_perm(grb.y.size());
            for (int i = 0; i < grb.y.size(); ++i)
              r_til_perm(i) = grb.y(i) / std::sqrt(std::max(dinv(i), 1e-12));
            grb.y = lmm_c.Q * (Xb0_til + r_til_perm);
          }
        }

        GenePrepLm lm_b;
        GenePrepLmm lmm_b;
        GenePrepGlm glm_b;
        GenePrepGlmm glmm_b;
        prep_null(model, opt.fast, grb, &lm_b, &lmm_b, &glm_b, &glmm_b);

        double minp = 1.0;
        if (cache_gtil) {
          for (size_t i = 0; i < cached_gtil.size(); ++i) {
            const AssocHit h = test_lmm_gtil(lmm_b, cached_gtil[i], cached_maf[i], ws_b);
            if (std::isfinite(h.p) && h.p < minp) minp = h.p;
          }
        } else {
          for (const auto& gd : cached_dosage) {
            const double p = test_one_p(model, opt.fast, grb, gd, &lm_b, &lmm_b, &glm_b, &glmm_b);
            if (std::isfinite(p) && p < minp) minp = p;
          }
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

// Per-thread cis parallel over a genotype source G (defined in scan_cis.cpp,
// instantiated for PlinkBed and VcfSession). Each thread opens its own G and
// streams region queries for its chunk of genes.
template <typename G>
void run_cis_parallel(const Options& opt, PhenoData& ph, const CovData& cov,
                      const std::unordered_map<std::string, GeneLoc>& annot,
                      Eigen::MatrixXd* Kptr, bool need_k, bool need_lmm_basis,
                      Model model, ScopeOut& so, double pthr,
                      std::vector<GeneSummary>& summaries, const std::string& src);

} // namespace eqtl
