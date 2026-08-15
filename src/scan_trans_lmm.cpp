#include "eqtl/scan_trans.hpp"
#include "eqtl/plink_bed.hpp"
#include "eqtl/vcf_session.hpp"
#include <omp.h>

namespace eqtl {

// ---------------------------------------------------------------------------
// LMM-specific helpers
// ---------------------------------------------------------------------------
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
  if (g_wss < 1e-12) { h.p = 1.0; return h; }
  if (ws.XtDX.rows() != p1 || ws.XtDX.cols() != p1) ws.XtDX.resize(p1, p1);
  if (ws.XtDy.size() != p1) ws.XtDy.resize(p1);
  ws.XtDX.noalias() = ws.Xg.transpose() * dinv.asDiagonal() * ws.Xg;
  ws.XtDy.noalias() = ws.Xg.transpose() * (dinv.asDiagonal() * prep.y_til);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(ws.XtDX);
  if (ldlt.info() != Eigen::Success) { h.p = 1.0; return h; }
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
  h.r2 = (prep.rss_null > 1e-15) ? std::max(0.0, 1.0 - q / prep.rss_null) : 0.0;
  return h;
}

static void free_lmm_gene_raw(GeneLmmJob& j) {
  j.gr.K.resize(0, 0);
  j.gr.basis.Q.resize(0, 0);
  j.gr.basis.lambda.resize(0);
  j.gr.has_basis = false;
}

// ---------------------------------------------------------------------------
// SNP-outer LMM for trans/gw
// ---------------------------------------------------------------------------
template <typename G>
void scan_lmm_snp_outer(const Options& opt, G& geno, const MissPolicy& mp, double maf,
                        const std::string& scope, ScopeOut& out, double pthr,
                        std::vector<GeneLmmJob>& jobs) {
  if (jobs.empty()) return;

  bool same_keep = true;
  for (size_t i = 1; i < jobs.size(); ++i) {
    if (jobs[i].gr.keep != jobs[0].gr.keep) { same_keep = false; break; }
  }

  LmmBasis shared_basis;
  bool have_shared_basis = false;
  if (same_keep && !jobs.empty()) {
    if (jobs[0].gr.has_basis) {
      shared_basis = jobs[0].gr.basis;
      have_shared_basis = true;
    } else {
      const Eigen::MatrixXd& K0 = jobs[0].gr.K_ref ? *jobs[0].gr.K_ref : jobs[0].gr.K;
      if (K0.rows() > 0) {
        Eigen::MatrixXd K_use = K0;
        if (opt.fast) sparsify_grm(K_use, 1e-4);
        shared_basis = make_lmm_basis(K_use);
        have_shared_basis = true;
      }
    }
  }

  for (size_t ji = 0; ji < jobs.size(); ++ji) {
    auto& j = jobs[ji];
    if (have_shared_basis) j.prep = prep_lmm(j.gr.y, j.gr.X, shared_basis, opt.fast);
    else if (j.gr.has_basis) j.prep = prep_lmm(j.gr.y, j.gr.X, j.gr.basis, opt.fast);
    else j.prep = prep_lmm(j.gr.y, j.gr.X, j.gr.K_ref ? *j.gr.K_ref : j.gr.K, opt.fast);
    // jobs[0].prep.Q is the shared Q used at test time; drop per-gene Q copies for the rest
    if (same_keep && ji > 0) j.prep.Q.resize(0, 0);
    if (opt.perm == 0) free_lmm_gene_raw(j);
    else {
      j.gr.K.resize(0, 0);
      // stage-2 uses the shared basis; drop the per-gene copy to avoid O(G·n²) memory
      if (have_shared_basis) {
        j.gr.basis.Q.resize(0, 0);
        j.gr.basis.lambda.resize(0);
        j.gr.has_basis = false;
      }
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

    size_t snp_cnt = 0;
    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      if (++snp_cnt % 100000 == 0)
        info("trans/gw LMM: SNP " + std::to_string(snp_cnt));
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
#pragma omp parallel if (opt.threads > 1 && jobs.size() > 32 && !omp_in_parallel()) num_threads(opt.threads)
      {
        LmmTestWs ws_t;
#pragma omp for schedule(static)
        for (int ji = 0; ji < Gz; ++ji) {
          auto& job = jobs[static_cast<size_t>(ji)];
          if (scope == "trans" && job.has_loc && in_cis_window(snp, job.loc, opt.window)) continue;
          AssocHit h = test_lmm_gtil(job.prep, g_til, maf_sub, ws_t);
          if (topK > 0 && std::isfinite(h.p) && h.p < opt.perm_trans_thr) topk_consider(job.top, topK, h.p, g_buf);
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
    size_t maxk = 0;
    for (const auto& j : jobs) maxk = std::max(maxk, j.gr.keep.size());
    g_buf.resize(static_cast<int>(maxk));
    g_til.resize(static_cast<int>(maxk));
    size_t snp_cnt2 = 0;
    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      if (++snp_cnt2 % 100000 == 0)
        info("trans/gw LMM: SNP " + std::to_string(snp_cnt2));
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
        {
          Eigen::Map<const Eigen::VectorXd> g(g_buf.data(), nk);
          g_til.noalias() = job.prep.Q.transpose() * g;
        }
        AssocHit h = test_lmm_gtil(job.prep, g_til, maf_sub, ws);
        if (opt.perm > 0 && std::isfinite(h.p) && h.p < opt.perm_trans_thr)
          topk_consider(job.top, opt.perm_trans_top, h.p,
                        Eigen::Map<const Eigen::VectorXd>(g_buf.data(), nk));
        apply_snp_hit(job, h, snp, pthr, Model::Lmm, out);
      }
      return true;
    });
  }

  for (auto& job : jobs) {
    job.summary.acat_p = job.acat_acc.p();
    if (job.best.p <= pthr && job.best.p <= 1.0)
      write_pair_line(out.top, job.best, Model::Lmm, scope);
    if (opt.perm > 0) {
      // Mixed keeps: each job's stage-2 needs its own basis (gr.K was freed above
      // and gr.basis never set, so prep_null would hit an empty K and abort).
      // prep.Q/lambda hold this job's decomposition and are dead after the scan:
      // move them into gr.basis so permutations reuse the basis instead of
      // re-decomposing per permutation.
      if (!have_shared_basis && !job.gr.has_basis) {
        job.gr.basis.Q = std::move(job.prep.Q);
        job.gr.basis.lambda = std::move(job.prep.lambda);
        job.gr.has_basis = true;
      }
      stage2_perm_topk(opt, Model::Lmm, job, have_shared_basis ? &shared_basis : nullptr);
    } else {
      job.summary.p_emp = std::numeric_limits<double>::quiet_NaN();
      job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

// Explicit instantiations
template void scan_lmm_snp_outer<PlinkBed>(const Options&, PlinkBed&, const MissPolicy&, double,
                                           const std::string&, ScopeOut&, double, std::vector<GeneLmmJob>&);
template void scan_lmm_snp_outer<VcfSession>(const Options&, VcfSession&, const MissPolicy&, double,
                                             const std::string&, ScopeOut&, double, std::vector<GeneLmmJob>&);

} // namespace eqtl
