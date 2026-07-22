#include "eqtl/scan.hpp"
#include "eqtl/output.hpp"
#include "eqtl/stats_extra.hpp"
#include "eqtl/util.hpp"
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

AssocHit run_test(Model model, bool fast, const Eigen::VectorXd& y, const Eigen::MatrixXd& X,
                  const Eigen::MatrixXd* K, const Eigen::VectorXd& g, GenePrepLm* lm_cache,
                  GenePrepLmm* lmm_cache, GenePrepGlm* glm_cache, GenePrepGlmm* glmm_cache,
                  bool have_cache) {
  switch (model) {
    case Model::Lm: {
      if (!have_cache) {
        *lm_cache = prep_lm(y, X);
      }
      return test_lm(*lm_cache, g);
    }
    case Model::Lmm: {
      if (!have_cache) {
        *lmm_cache = prep_lmm(y, X, *K, fast);
      }
      return test_lmm(*lmm_cache, g);
    }
    case Model::Glm: {
      if (!have_cache) {
        *glm_cache = prep_glm_nb(y, X, fast);
      }
      return test_glm_nb(*glm_cache, g);
    }
    case Model::Glmm: {
      if (!have_cache) {
        *glmm_cache = prep_glmm_pois(y, X, *K, fast);
      }
      return test_glmm_pois(*glmm_cache, g);
    }
  }
  return {};
}

void scan_gene_snps(const Options& opt, Model model, const std::string& scope, const std::string& gene,
                    const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd* K,
                    const GeneLoc* loc, const std::vector<SnpRec>& snps, double pthr, ScopeOut& out,
                    GeneSummary& summary) {
  GenePrepLm lm_c;
  GenePrepLmm lmm_c;
  GenePrepGlm glm_c;
  GenePrepGlmm glmm_c;
  bool cache = false;

  AssocHit best;
  best.p = 2.0;
  std::vector<double> pvals;
  pvals.reserve(snps.size());

  // prebuild gene prep once
  Eigen::VectorXd g0 = Eigen::VectorXd::Zero(y.size());
  run_test(model, opt.fast, y, X, K, g0, &lm_c, &lmm_c, &glm_c, &glmm_c, false);
  cache = true;

  for (const auto& snp : snps) {
    Eigen::Map<const Eigen::VectorXd> g(snp.dosage.data(), static_cast<int>(snp.dosage.size()));
    AssocHit h = run_test(model, opt.fast, y, X, K, g, &lm_c, &lmm_c, &glm_c, &glmm_c, cache);
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
    if (std::isfinite(h.p)) {
      pvals.push_back(h.p);
    }
    if (h.p <= pthr) {
      write_pair_line(out.pairs, h, model, scope);
      ++summary.n_sig;
    }
    if (h.p < best.p) {
      best = h;
    }
  }

  summary.gene = gene;
  summary.n_tested = static_cast<int>(snps.size());
  if (loc && loc->ok) {
    summary.chrom = loc->chrom;
    summary.tss = loc->tss;
  } else if (!snps.empty()) {
    summary.chrom = snps.front().chrom;
  }
  summary.acat_p = acat(pvals);

  // empirical p (tensorqtl-class)
  if (opt.perm > 0 && !snps.empty()) {
    std::mt19937 rng(opt.seed >= 0 ? opt.seed : 1);
    // phenotype permutation: shuffle y indices
    std::vector<int> idx(y.size());
    std::iota(idx.begin(), idx.end(), 0);
    const double T_obs = best.p;  // min p
    std::vector<double> T_perm;
    T_perm.reserve(opt.perm);
    std::vector<double> perm_min_p;
    perm_min_p.reserve(opt.perm);

    for (int b = 0; b < opt.perm; ++b) {
      std::shuffle(idx.begin(), idx.end(), rng);
      Eigen::VectorXd yb(y.size());
      for (int i = 0; i < y.size(); ++i) {
        yb(i) = y(idx[i]);
      }
      GenePrepLm lm_b;
      GenePrepLmm lmm_b;
      GenePrepGlm glm_b;
      GenePrepGlmm glmm_b;
      Eigen::VectorXd g0b = Eigen::VectorXd::Zero(y.size());
      run_test(model, opt.fast, yb, X, K, g0b, &lm_b, &lmm_b, &glm_b, &glmm_b, false);
      double minp = 1.0;
      for (const auto& snp : snps) {
        Eigen::Map<const Eigen::VectorXd> g(snp.dosage.data(), static_cast<int>(snp.dosage.size()));
        AssocHit hb = run_test(model, opt.fast, yb, X, K, g, &lm_b, &lmm_b, &glm_b, &glmm_b, true);
        if (hb.p < minp) {
          minp = hb.p;
        }
      }
      // T larger is more extreme: use -log10(p)
      T_perm.push_back(-std::log10(std::max(minp, 1e-300)));
      perm_min_p.push_back(minp);
    }
    const double Tobs = -std::log10(std::max(T_obs, 1e-300));
    summary.p_emp = p_emp_count(Tobs, T_perm);
    if (!opt.disable_beta_approx) {
      beta_approx_p(perm_min_p, T_obs, summary.p_beta, summary.beta_shape1, summary.beta_shape2);
    } else {
      summary.p_beta = summary.p_emp;
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

std::vector<SnpRec> collect_cis(VcfSession& vcf, const GeneLoc& loc, int window, MissHand miss) {
  std::vector<SnpRec> snps;
  if (!loc.ok) {
    return snps;
  }
  const int64_t start = std::max<int64_t>(1, loc.tss - window);
  const int64_t end = loc.tss + window;
  vcf.for_each_snp_region(loc.chrom, start, end, miss, [&](const SnpRec& s) {
    snps.push_back(s);
    return true;
  });
  return snps;
}

}  // namespace

int run_make_grm(const Options& opt) {
  VcfSession vcf;
  vcf.open(opt.vcf);
  std::vector<std::string> samples = vcf.samples();
  if (!opt.pheno.empty()) {
    PhenoData ph = load_pheno(opt.pheno);
    samples = intersect_order(ph.sample_ids, vcf.samples());
    if (samples.empty()) {
      die("no overlapping samples for GRM");
    }
  }
  Grm g = compute_grm(vcf, samples, opt.miss);
  write_grm_gcta(opt.out, g);
  return 0;
}

int run_eqtl(const Options& opt) {
  omp_set_num_threads(opt.threads);

  PhenoData ph = load_pheno(opt.pheno);
  VcfSession vcf;
  vcf.open(opt.vcf);

  std::vector<std::string> sample_order = intersect_order(ph.sample_ids, vcf.samples());
  if (sample_order.empty()) {
    die("no overlapping samples between pheno and VCF");
  }
  info("overlap samples: " + std::to_string(sample_order.size()));
  slice_pheno(ph, sample_order);
  CovData cov = load_covar(opt.covar, sample_order);
  vcf.set_sample_order(sample_order);

  std::unordered_map<std::string, GeneLoc> annot;
  const bool have_gff = !opt.gff.empty();
  if (have_gff) {
    annot = load_gff_tss(opt.gff, opt.gff_id_key);
  } else {
    info("no GFF: mode ignored; genome-wide all-pairs (gw)");
  }

  // determine effective modes
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
    if (needs_grm(m)) {
      any_grm = true;
    }
  }
  Grm grm;
  Eigen::MatrixXd* Kptr = nullptr;
  Eigen::MatrixXd Kmat;
  if (any_grm) {
    if (!opt.grm.empty()) {
      grm = slice_grm(load_grm_gcta(opt.grm), sample_order);
    } else {
      info("computing GRM from VCF (overlap samples)");
      grm = compute_grm(vcf, sample_order, opt.miss);
    }
    Kmat = grm.K;
    Kptr = &Kmat;
  }

  // For smoke/v1: load SNPs into memory per scope strategy
  // IO-max: stream full VCF once for gw/trans; cis uses region queries
  std::vector<SnpRec> all_snps;
  const bool need_all = std::find(scopes.begin(), scopes.end(), "trans") != scopes.end() ||
                        std::find(scopes.begin(), scopes.end(), "gw") != scopes.end();
  if (need_all) {
    info("loading SNPs (stream)");
    all_snps = vcf.load_all(opt.miss, -1);
    info("SNPs loaded: " + std::to_string(all_snps.size()));
  }

  for (Model model : opt.models) {
    if (needs_grm(model) && !Kptr) {
      die("internal: GRM required");
    }
    for (const std::string& scope : scopes) {
      const std::string prefix = opt.out + "." + model_str(model) + "." + scope;
      ScopeOut so;
      so.tag = scope;
      so.pairs.open(prefix + ".pairs.tsv");
      so.top.open(prefix + ".top.tsv");
      so.region.open(prefix + ".region.tsv");
      if (!so.pairs || !so.top || !so.region) {
        die("cannot open output for " + prefix);
      }
      write_pairs_header(so.pairs, model);
      write_top_header(so.top, model);
      write_region_header(so.region);

      const double pthr = (scope == "cis") ? opt.pval_cis : opt.pval_trans;

      for (size_t gi = 0; gi < ph.gene_ids.size(); ++gi) {
        const std::string& gene = ph.gene_ids[gi];
        Eigen::VectorXd y = ph.Y.col(static_cast<int>(gi));

        // skip genes with NA / zero variance
        if (!y.allFinite()) {
          continue;
        }
        if (y.array().abs().maxCoeff() < 1e-15) {
          continue;
        }
        if (needs_counts(model) && !looks_like_counts(y)) {
          warn("skip non-count gene for " + model_str(model) + ": " + gene);
          continue;
        }

        const GeneLoc* locp = nullptr;
        GeneLoc loc_store;
        if (have_gff) {
          auto it = annot.find(gene);
          if (it == annot.end()) {
            continue;  // 0 hit skip
          }
          loc_store = it->second;
          locp = &loc_store;
        }

        std::vector<SnpRec> snps;
        if (scope == "cis") {
          if (!locp) {
            continue;
          }
          snps = collect_cis(vcf, *locp, opt.window, opt.miss);
        } else if (scope == "trans") {
          if (!locp) {
            continue;
          }
          const int64_t cstart = std::max<int64_t>(1, locp->tss - opt.window);
          const int64_t cend = locp->tss + opt.window;
          for (const auto& s : all_snps) {
            if (s.chrom == locp->chrom && s.pos >= cstart && s.pos <= cend) {
              continue;
            }
            snps.push_back(s);
          }
        } else {  // gw
          snps = all_snps;
        }

        if (snps.empty()) {
          continue;
        }

        GeneSummary summary;
        scan_gene_snps(opt, model, scope, gene, y, cov.X, Kptr, locp, snps, pthr, so, summary);
      }
      info("finished " + prefix);
    }
  }
  return 0;
}

}  // namespace eqtl
