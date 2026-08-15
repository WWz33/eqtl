/* eqtl — scan orchestrator: dispatch cis/trans/gw × model × genotype source */
#include "eqtl/scan.hpp"
#include "eqtl/scan_common.hpp"
#include "eqtl/scan_cis.hpp"
#include "eqtl/scan_trans.hpp"
#include "eqtl/plink_bed.hpp"
#include "eqtl/vcf_session.hpp"
#include <omp.h>
#include <cstdio>
#include <unordered_set>
#include <unordered_map>

namespace eqtl {

// ---------------------------------------------------------------------------
// --make-grm
// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// Main scan driver (templated on genotype source)
// ---------------------------------------------------------------------------
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

  // --- upfront validation (fail fast before expensive GRM / scan) ---
  {
    std::ofstream probe(opt.out + ".probe");
    if (!probe) die("cannot write to output prefix: " + opt.out);
    probe.close();
    std::remove((opt.out + ".probe").c_str());
  }
  {
    bool need_region = false;
    for (const auto& s : scopes)
      if (s == "cis") need_region = true;
    if (need_region && !geno.has_index())
      die("cis mode requires an indexed genotype file; for VCF run: tabix -p vcf <file>.vcf.gz");
  }
  if (have_gff) {
    size_t matched = 0;
    for (const auto& g : ph.gene_ids)
      if (annot.count(g)) ++matched;
    info("pheno genes matched to GFF: " + std::to_string(matched) + "/" +
         std::to_string(ph.gene_ids.size()));
    if (matched == 0)
      die("no phenotype gene IDs match the GFF annotation (check --gff-id-key and ID format)");
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
            if (it != annot.end()) { loc_store = it->second; has_loc = true; }
            else if (scope != "gw") continue;
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
            if (it != annot.end()) { loc_store = it->second; has_loc = true; }
            else if (scope != "gw") continue;
          }
          if (scope == "trans" && !has_loc) continue;
          GeneLmmJob job;
          job.gene = gene;
          job.loc = loc_store;
          job.has_loc = has_loc;
          if (!build_gene_ready(y, cov.X, Kptr, need_k, false, opt.fast, job.gr)) continue;
          jobs.push_back(std::move(job));
        }
        info("trans/gw LMM: SNP-outer (" + std::to_string(jobs.size()) + " genes)");
        scan_lmm_snp_outer(opt, geno, mp, maf, scope, so, pthr, jobs);
        summaries.reserve(jobs.size());
        for (auto& j : jobs) summaries.push_back(std::move(j.summary));

      } else if (scope == "cis" && opt.threads > 1 && have_gff) {
        info("cis: per-thread parallel (" + std::to_string(opt.threads) + " threads)");
        if (opt.use_bfile()) {
          run_cis_parallel<PlinkBed>(opt, ph, cov, annot, Kptr, need_k, need_lmm_basis, model, so,
                                     pthr, summaries, opt.bfile);
        } else {
          run_cis_parallel<VcfSession>(opt, ph, cov, annot, Kptr, need_k, need_lmm_basis, model, so,
                                       pthr, summaries, opt.vcf);
        }

      } else {
        // Gene-outer fallback (cis single-thread, GLM/GLMM, trans/gw non-LM/LMM)
        LmmBasis grm_basis;
        bool have_grm_basis = false;
        if (need_lmm_basis && Kptr && Kptr->rows() == static_cast<int>(ph.sample_ids.size())) {
          Eigen::MatrixXd K_use = *Kptr;
          if (opt.fast) sparsify_grm(K_use, 1e-4);
          grm_basis = make_lmm_basis(K_use);
          have_grm_basis = true;
          info("LMM: shared GRM eigen-decomp (n=" + std::to_string(Kptr->rows()) + ")");
        }

        for (size_t gi = 0; gi < ph.gene_ids.size(); ++gi) {
          if (gi % 500 == 0)
            info("scan: gene " + std::to_string(gi) + "/" + std::to_string(ph.gene_ids.size()));
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
            if (it != annot.end()) { loc_store = it->second; locp = &loc_store; }
            else if (scope != "gw") continue;
          }

          GeneReady gr;
          if (have_grm_basis && all_finite(y, cov.X)) {
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
            const auto [cstart, cend] = cis_window_bounds(*locp, opt.window);
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

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------
int run_eqtl(const Options& opt) {
  omp_set_num_threads(opt.threads);
  omp_set_max_active_levels(2);

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

} // namespace eqtl
