#include "eqtl/scan.hpp"
#include "eqtl/output.hpp"
#include "eqtl/stats_extra.hpp"
#include "eqtl/util.hpp"
#include "eqtl/plink_bed.hpp"
#include <fstream>
#include <random>
#include <omp.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_map>

namespace eqtl {

namespace {

struct ScopeOut {
  std::ofstream pairs;
  std::ofstream top;
  std::ofstream region;
  std::string tag;
};

bool needs_grm(Model m) {
  return m == Model::Lmm || m == Model::Glmm;
}

bool needs_counts(Model m) {
  return m == Model::Glm || m == Model::Glmm;
}

// GCTA-style: keep individuals with finite y (and finite covar rows)
struct GeneReady {
  std::vector<int> keep; // indices into full sample order
  Eigen::VectorXd y;
  Eigen::MatrixXd X;
  Eigen::MatrixXd K; // n_keep x n_keep if used
  LmmBasis basis;
  bool has_basis = false;
  int n_full = 0;
};

bool build_gene_ready(const Eigen::VectorXd& y_full, const Eigen::MatrixXd& X_full,
                      const Eigen::MatrixXd* K_full, bool need_k, bool need_lmm_basis,
                      bool fast_sparse, GeneReady& out) {
  out.n_full = static_cast<int>(y_full.size());
  out.keep.clear();
  out.keep.reserve(out.n_full);
  for (int i = 0; i < out.n_full; ++i) {
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
      Eigen::MatrixXd K_use = out.K;
      if (fast_sparse) sparsify_grm(K_use, 1e-4);
      out.basis = make_lmm_basis(K_use);
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

// If dosage length already equals keep size (pre-subset), use as-is
Eigen::VectorXd dosage_vec(const SnpRec& snp, const GeneReady& gr) {
  if (static_cast<int>(snp.dosage.size()) == static_cast<int>(gr.keep.size())) {
    return Eigen::Map<const Eigen::VectorXd>(snp.dosage.data(), static_cast<int>(snp.dosage.size()));
  }
  return subset_dosage(snp.dosage, gr.keep);
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

// Test one SNP; returns hit with meta filled from snp/gene
AssocHit test_one(Model model, bool fast, const GeneReady& gr, const SnpRec& snp,
                  const std::string& gene, const GeneLoc* loc, GenePrepLm* lm_c, GenePrepLmm* lmm_c,
                  GenePrepGlm* glm_c, GenePrepGlmm* glmm_c, bool have_cache) {
  Eigen::VectorXd g = dosage_vec(snp, gr);
  // drop SNP if any non-finite dosage among keep (after impute path should be finite)
  for (int i = 0; i < g.size(); ++i) {
    if (!std::isfinite(g(i))) {
      AssocHit h;
      h.p = 1.0;
      return h;
    }
  }
  AssocHit h = run_test(model, fast, gr, g, lm_c, lmm_c, glm_c, glmm_c, have_cache);
  h.gene = gene;
  h.snp = snp.id;
  h.chrom = snp.chrom;
  h.pos = snp.pos;
  h.ref = snp.ref;
  h.alt = snp.alt;
  h.maf = snp.maf;
  if (loc && loc->ok) {
    h.has_tss_dist = true;
    h.tss_dist = static_cast<double>(snp.pos - loc->tss);
  }
  return h;
}

bool in_cis_window(const SnpRec& s, const GeneLoc& loc, int window) {
  if (!chrom_equal(s.chrom, loc.chrom)) return false;
  const int64_t cstart = std::max<int64_t>(1, loc.tss - window);
  const int64_t cend = loc.tss + window;
  return s.pos >= cstart && s.pos <= cend;
}

// Collect SNP list only for small cis regions; null for stream modes
void scan_gene_snps(const Options& opt, Model model, const std::string& scope, const std::string& gene,
                    const GeneReady& gr, const GeneLoc* loc, const std::vector<SnpRec>* snps_list,
                    // if snps_list null and stream_fn set, caller streams; here we only handle list
                    double pthr, ScopeOut& out, GeneSummary& summary,
                    const std::function<void(const std::function<void(const SnpRec&)>&)>& stream_snps) {
  GenePrepLm lm_c;
  GenePrepLmm lmm_c;
  GenePrepGlm glm_c;
  GenePrepGlmm glmm_c;

  // build null cache once
  {
    Eigen::VectorXd g0 = Eigen::VectorXd::Zero(gr.y.size());
    run_test(model, opt.fast, gr, g0, &lm_c, &lmm_c, &glm_c, &glmm_c, false);
  }

  AssocHit best;
  best.p = 2.0;
  std::vector<double> pvals;
  int n_tested = 0;
  summary.n_sig = 0;

  auto consume = [&](const SnpRec& snp) {
    AssocHit h =
        test_one(model, opt.fast, gr, snp, gene, loc, &lm_c, &lmm_c, &glm_c, &glmm_c, true);
    if (!std::isfinite(h.p)) return;
    ++n_tested;
    pvals.push_back(h.p);
    if (h.p <= pthr) {
      write_pair_line(out.pairs, h, model, scope);
      ++summary.n_sig;
    }
    if (h.p < best.p) best = h;
  };

  if (snps_list) {
    const int n = static_cast<int>(snps_list->size());
    std::vector<AssocHit> hits(static_cast<size_t>(n));
#pragma omp parallel for schedule(static) if (opt.threads > 1)
    for (int si = 0; si < n; ++si) {
      hits[static_cast<size_t>(si)] = test_one(model, opt.fast, gr, (*snps_list)[static_cast<size_t>(si)],
                                               gene, loc, &lm_c, &lmm_c, &glm_c, &glmm_c, true);
    }
    for (const auto& h : hits) {
      if (!std::isfinite(h.p)) continue;
      ++n_tested;
      pvals.push_back(h.p);
      if (h.p <= pthr) {
        write_pair_line(out.pairs, h, model, scope);
        ++summary.n_sig;
      }
      if (h.p < best.p) best = h;
    }
  } else if (stream_snps) {
    stream_snps(consume);
  }

  summary.gene = gene;
  summary.n_tested = n_tested;
  if (loc && loc->ok) {
    summary.chrom = loc->chrom;
    summary.tss = loc->tss;
  }
  summary.acat_p = acat(pvals);

  // gene-level perm: shuffle y among keep; re-fit null; min p
  if (opt.perm > 0 && n_tested > 0) {
    const double T_obs = best.p;
    std::vector<double> T_perm(static_cast<size_t>(opt.perm));
    std::vector<double> perm_min_p(static_cast<size_t>(opt.perm));

    const bool have_list = snps_list != nullptr;

#pragma omp parallel for schedule(dynamic) if (opt.threads > 1)
    for (int b = 0; b < opt.perm; ++b) {
      std::mt19937 rng_b(static_cast<unsigned>(opt.seed >= 0 ? opt.seed : 1) +
                         static_cast<unsigned>(b) * 9973u +
                         static_cast<unsigned>(std::hash<std::string>{}(gene)));

      std::vector<int> idx(static_cast<size_t>(gr.y.size()));
      std::iota(idx.begin(), idx.end(), 0);
      std::shuffle(idx.begin(), idx.end(), rng_b);

      GeneReady grb = gr;
      for (int i = 0; i < gr.y.size(); ++i) grb.y(i) = gr.y(idx[static_cast<size_t>(i)]);

      GenePrepLm lm_b;
      GenePrepLmm lmm_b;
      GenePrepGlm glm_b;
      GenePrepGlmm glmm_b;
      Eigen::VectorXd g0 = Eigen::VectorXd::Zero(grb.y.size());
      run_test(model, opt.fast, grb, g0, &lm_b, &lmm_b, &glm_b, &glmm_b, false);

      double minp = 1.0;
      auto take = [&](const SnpRec& snp) {
        AssocHit hb =
            test_one(model, opt.fast, grb, snp, gene, loc, &lm_b, &lmm_b, &glm_b, &glmm_b, true);
        if (std::isfinite(hb.p) && hb.p < minp) minp = hb.p;
      };

      if (have_list) {
        const int n = static_cast<int>(snps_list->size());
        for (int si = 0; si < n; ++si) take((*snps_list)[static_cast<size_t>(si)]);
      } else if (stream_snps) {
#pragma omp critical(eqtl_perm_stream)
        {
          stream_snps(take);
        }
      }

      T_perm[static_cast<size_t>(b)] = -std::log10(std::max(minp, 1e-300));
      perm_min_p[static_cast<size_t>(b)] = minp;
    }

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
  write_region_line(out.region, summary);
}

}  // namespace

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

template <typename G>
int run_eqtl_geno(const Options& opt, G& geno, PhenoData& ph,
                  const std::vector<std::string>& sample_order) {
  CovData cov = load_covar(opt.covar, sample_order);
  geno.set_sample_order(sample_order);

  std::unordered_map<std::string, GeneLoc> annot;
  const bool have_gff = !opt.gff.empty();
  if (have_gff) {
    annot = load_gff_tss(opt.gff, opt.gff_id_key);
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
  info(have_gff && scopes.size() == 1 && scopes[0] == "cis" ? "cis-only: region queries"
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
      static char pairs_buf[1 << 20];
      static char top_buf[1 << 16];
      static char region_buf[1 << 16];
      so.pairs.rdbuf()->pubsetbuf(pairs_buf, sizeof(pairs_buf));
      so.top.rdbuf()->pubsetbuf(top_buf, sizeof(top_buf));
      so.region.rdbuf()->pubsetbuf(region_buf, sizeof(region_buf));
      write_pairs_header(so.pairs, model);
      write_top_header(so.top, model);
      write_region_header(so.region);

      const double pthr = (scope == "cis") ? opt.pval_cis : opt.pval_trans;

      for (size_t gi = 0; gi < ph.gene_ids.size(); ++gi) {
        const std::string& gene = ph.gene_ids[gi];
        Eigen::VectorXd y = ph.Y.col(static_cast<int>(gi));

        if (needs_counts(model) && !looks_like_counts(y)) {
          // allow NA: check on finite subset later
          bool any_count = false;
          int checked = 0;
          for (int i = 0; i < y.size() && checked < 20; ++i) {
            if (!std::isfinite(y(i))) continue;
            ++checked;
            if (y(i) < 0 || std::fabs(y(i) - std::floor(y(i))) > 1e-8) {
              any_count = false;
              break;
            }
            any_count = true;
          }
          if (!any_count) {
            warn("skip non-count gene for " + model_str(model) + ": " + gene);
            continue;
          }
        }

        const GeneLoc* locp = nullptr;
        GeneLoc loc_store;
        if (have_gff) {
          auto it = annot.find(gene);
          if (it == annot.end()) continue;
          loc_store = it->second;
          locp = &loc_store;
        }

        GeneReady gr;
        if (!build_gene_ready(y, cov.X, Kptr, need_k, need_lmm_basis, opt.fast, gr)) continue;

        GeneSummary summary;

        if (scope == "cis") {
          if (!locp) continue;
          const int64_t cstart = std::max<int64_t>(1, locp->tss - opt.window);
          const int64_t cend = locp->tss + opt.window;
          std::vector<SnpRec> cis_snps = geno.load_region(locp->chrom, cstart, cend, mp, maf);
          scan_gene_snps(opt, model, scope, gene, gr, locp, &cis_snps, pthr, so, summary, nullptr);
        } else if (scope == "trans") {
          if (!locp) continue;
          auto stream = [&](const std::function<void(const SnpRec&)>& take) {
            geno.for_each_snp(mp, maf, [&](const SnpRec& s) {
              if (in_cis_window(s, *locp, opt.window)) return true;
              take(s);
              return true;
            });
          };
          scan_gene_snps(opt, model, scope, gene, gr, locp, nullptr, pthr, so, summary, stream);
        } else {  // gw
          auto stream = [&](const std::function<void(const SnpRec&)>& take) {
            geno.for_each_snp(mp, maf, [&](const SnpRec& s) {
              take(s);
              return true;
            });
          };
          scan_gene_snps(opt, model, scope, gene, gr, locp, nullptr, pthr, so, summary, stream);
        }
      }
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
  std::vector<std::string> sample_order = intersect_order(ph.sample_ids, vcf.samples());
  if (sample_order.empty()) die("no overlapping samples between pheno and VCF");
  info("overlap samples: " + std::to_string(sample_order.size()));
  slice_pheno(ph, sample_order);
  return run_eqtl_geno(opt, vcf, ph, sample_order);
}

}  // namespace eqtl
