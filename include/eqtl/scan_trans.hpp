/* eqtl — SNP-outer trans/gw scan: job types, helpers, template declarations */
#pragma once
#include "eqtl/scan_common.hpp"
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>

namespace eqtl {

// ---------------------------------------------------------------------------
// Job types for SNP-outer trans/gw
// ---------------------------------------------------------------------------
struct GeneLmJob {
  std::string gene;
  GeneLoc loc;
  bool has_loc = false;
  GeneReady gr;
  GenePrepLm prep;
  GeneSummary summary;
  AssocHit best;
  AcatAcc acat_acc;
  std::vector<std::pair<double, Eigen::VectorXd>> top;
};

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

struct LmmTestWs {
  Eigen::MatrixXd Xg;
  Eigen::MatrixXd XtDX;
  Eigen::VectorXd XtDy, beta, e, cov_col;
};

// ---------------------------------------------------------------------------
// Inline helpers shared by LM/LMM SNP-outer paths
// ---------------------------------------------------------------------------

// max-heap by p (largest p at front) — retain smallest K p-values
inline void topk_consider(std::vector<std::pair<double, Eigen::VectorXd>>& top, int K, double p,
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
inline void fill_hit_meta(AssocHit& h, const Job& job, const SnpRec& snp) {
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

template <typename Job>
inline bool apply_snp_hit_stats(Job& job, AssocHit& h, const SnpRec& snp, double pthr) {
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
inline void apply_snp_hit(Job& job, AssocHit& h, const SnpRec& snp, double pthr, Model model,
                          ScopeOut& out) {
  if (apply_snp_hit_stats(job, h, snp, pthr)) {
    write_pair_line(out.pairs, h, model, out.tag);
  }
}

// ---------------------------------------------------------------------------
// Template declarations (defined in scan_trans_lm.cpp / scan_trans_lmm.cpp)
// ---------------------------------------------------------------------------
template <typename G>
void scan_lm_snp_outer(const Options& opt, G& geno, const MissPolicy& mp, double maf,
                       const std::string& scope, ScopeOut& out, double pthr,
                       std::vector<GeneLmJob>& jobs);

template <typename G>
void scan_lmm_snp_outer(const Options& opt, G& geno, const MissPolicy& mp, double maf,
                        const std::string& scope, ScopeOut& out, double pthr,
                        std::vector<GeneLmmJob>& jobs);

// Stage-2 gene perm on top-K SNPs (defined in scan_perm.cpp)
template <typename Job>
void stage2_perm_topk(const Options& opt, Model model, Job& job,
                      const LmmBasis* ext_basis = nullptr);

} // namespace eqtl
