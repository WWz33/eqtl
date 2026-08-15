#include "eqtl/scan_cis.hpp"
#include "eqtl/plink_bed.hpp"
#include "eqtl/vcf_session.hpp"
#include <omp.h>
#include <mutex>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <sstream>

namespace eqtl {

// Genes are already parallelized across threads; drop per-thread BGZF decode
// threads so each thread keeps a single decode context.
static void prep_thread_geno(PlinkBed&) {}
static void prep_thread_geno(VcfSession& v) { v.set_threads(1); }

template <typename G>
void run_cis_parallel(const Options& opt, PhenoData& ph, const CovData& cov,
                      const std::unordered_map<std::string, GeneLoc>& annot,
                      Eigen::MatrixXd* Kptr, bool need_k, bool need_lmm_basis,
                      Model model, ScopeOut& so, double pthr,
                      std::vector<GeneSummary>& summaries, const std::string& src) {
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
  std::string err_msg;
  std::mutex err_mu;

  LmmBasis grm_basis;
  bool have_grm_basis = false;
  if (need_lmm_basis && Kptr && Kptr->rows() == static_cast<int>(ph.sample_ids.size())) {
    Eigen::MatrixXd K_use = *Kptr;
    if (opt.fast) sparsify_grm(K_use, 1e-4);
    grm_basis = make_lmm_basis(K_use);
    have_grm_basis = true;
  }

  // Pre-open the genotype source once before the parallel region. For PlinkBed
  // this amortizes the BIM/FAM parse (~tens of MB of string tables) across all
  // worker threads via clone_for_thread(meta-only shared, per-thread bed FILE*).
  // For VcfSession each thread still opens its own htsFile/header — opening one
  // here only to share metadata is not possible for the indexed-BCF VCF path.
  std::unique_ptr<G> g_master;
  if constexpr (std::is_same_v<G, PlinkBed>) {
    g_master = std::make_unique<G>();
    g_master->open(src);
  }

#pragma omp parallel num_threads(T)
  {
    const int tid = omp_get_thread_num();
    const int nT = omp_get_num_threads();
    try {
      G g;
      if constexpr (std::is_same_v<G, PlinkBed>) {
        if (g_master) PlinkBed::clone_for_thread(*g_master, g);
        else g.open(src);
      } else {
        g.open(src);
      }
      prep_thread_geno(g);
      g.set_sample_order(ph.sample_ids);
      ScopeOut local;
      const std::string t_pairs = opt.out + ".tmp." + std::to_string(tid) + ".pairs";
      const std::string t_top = opt.out + ".tmp." + std::to_string(tid) + ".top";
      const std::string t_region = opt.out + ".tmp." + std::to_string(tid) + ".region";
      local.tag = "cis";
      local.pairs.open(t_pairs);
      local.top.open(t_top);
      local.region.open(t_region);

      const size_t nW = works.size();
      const size_t chunk = (nW + static_cast<size_t>(nT) - 1) / static_cast<size_t>(nT);
      const size_t lo = static_cast<size_t>(tid) * chunk;
      const size_t hi = std::min(nW, lo + chunk);
      for (size_t wi = lo; wi < hi; ++wi) {
        if (err.load()) break;
        if (tid == 0 && ((wi - lo) % 100 == 0))
          info("cis: gene " + std::to_string(wi) + "/" + std::to_string(nW));
        const auto& w = works[wi];
        Eigen::VectorXd y = ph.Y.col(w.gi);
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
        GeneSummary summary;
        const GeneLoc* locp = &w.loc;
        const auto [cstart, cend] = cis_window_bounds(*locp, opt.window);
        auto stream = [&](auto&& take) {
          g.for_each_snp_region(locp->chrom, cstart, cend, mp, maf, [&](const SnpRec& s) {
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
      (void)t_region;
    } catch (const std::exception& e) {
      std::lock_guard<std::mutex> lk(err_mu);
      if (err_msg.empty()) err_msg = e.what();
      err.store(1);
    } catch (...) {
      err.store(1);
    }
  }
  if (err.load())
    die(err_msg.empty() ? "cis parallel scan failed" : "cis parallel scan failed: " + err_msg);

  auto cat_file = [](std::ostream& dst, const std::string& path) {
    std::ifstream in(path);
    if (!in) return;
    dst << in.rdbuf();
    std::remove(path.c_str());
  };
  for (int t = 0; t < T; ++t) {
    cat_file(so.pairs, pairs_part[static_cast<size_t>(t)]);
    cat_file(so.top, top_part[static_cast<size_t>(t)]);
  }
  for (int t = 0; t < T; ++t)
    for (auto& s : part[static_cast<size_t>(t)])
      summaries.push_back(std::move(s));
}

template void run_cis_parallel<PlinkBed>(const Options&, PhenoData&, const CovData&,
                                         const std::unordered_map<std::string, GeneLoc>&,
                                         Eigen::MatrixXd*, bool, bool, Model, ScopeOut&, double,
                                         std::vector<GeneSummary>&, const std::string&);
template void run_cis_parallel<VcfSession>(const Options&, PhenoData&, const CovData&,
                                           const std::unordered_map<std::string, GeneLoc>&,
                                           Eigen::MatrixXd*, bool, bool, Model, ScopeOut&, double,
                                           std::vector<GeneSummary>&, const std::string&);

} // namespace eqtl
