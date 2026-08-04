#include "eqtl/scan_trans.hpp"
#include "eqtl/plink_bed.hpp"
#include "eqtl/vcf_session.hpp"
#include <omp.h>

namespace eqtl {

// ---------------------------------------------------------------------------
// LM-specific helpers
// ---------------------------------------------------------------------------
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
  if (gtg < 1e-12) { h.p = 1.0; return h; }
  h.beta = gty / gtg;
  const double df = prep.n - prep.p - 1;
  if (df <= 0) { h.p = 1.0; return h; }
  const double rss = prep.yty - h.beta * gty;
  const double s2 = std::max(rss / df, 0.0);
  h.se = std::sqrt(s2 / gtg);
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  h.r2 = (prep.yty > 0) ? (1.0 - rss / prep.yty) : 0.0;
  return h;
}

// ---------------------------------------------------------------------------
// SNP-outer LM for trans/gw
// ---------------------------------------------------------------------------
template <typename G>
void scan_lm_snp_outer(const Options& opt, G& geno, const MissPolicy& mp, double maf,
                       const std::string& scope, ScopeOut& out, double pthr,
                       std::vector<GeneLmJob>& jobs) {
  if (jobs.empty()) return;

  bool same_keep = true;
  for (size_t i = 1; i < jobs.size(); ++i) {
    if (jobs[i].gr.keep != jobs[0].gr.keep) { same_keep = false; break; }
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
    const int n = jobs[0].prep.n;
    const int Gz = static_cast<int>(jobs.size());
    Eigen::MatrixXd Ys(Gz, n);
    for (int gi = 0; gi < Gz; ++gi) Ys.row(gi) = jobs[static_cast<size_t>(gi)].prep.y_s.transpose();
    const auto& keep = jobs[0].gr.keep;
    const GenePrepLm& prep0 = jobs[0].prep;
    Eigen::VectorXd g_full(n), g_s(n), gty(Gz), Xt_g_buf(prep0.p);
    std::vector<AssocHit> write_hits(static_cast<size_t>(Gz));
    std::vector<char> write_flag(static_cast<size_t>(Gz), 0);
    const int topK = (opt.perm > 0) ? opt.perm_trans_top : 0;

    size_t snp_cnt = 0;
    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      if (++snp_cnt % 100000 == 0)
        info("trans/gw LM: SNP " + std::to_string(snp_cnt));
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
          h.p = 1.0; h.n = prep0.n; h.maf = maf_sub;
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
        if (topK > 0 && std::isfinite(h.p) && h.p < opt.perm_trans_thr) topk_consider(job.top, topK, h.p, g_full);
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
    size_t maxk = 0;
    for (const auto& j : jobs) maxk = std::max(maxk, j.gr.keep.size());
    Eigen::VectorXd g_buf(static_cast<int>(maxk));
    size_t snp_cnt2 = 0;
    geno.for_each_snp(mp, maf, [&](const SnpRec& snp) {
      if (++snp_cnt2 % 100000 == 0)
        info("trans/gw LM: SNP " + std::to_string(snp_cnt2));
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
        if (opt.perm > 0 && std::isfinite(h.p) && h.p < opt.perm_trans_thr)
          topk_consider(job.top, opt.perm_trans_top, h.p,
                        Eigen::Map<const Eigen::VectorXd>(g_buf.data(), nk));
        apply_snp_hit(job, h, snp, pthr, Model::Lm, out);
      }
      return true;
    });
  }

  for (auto& job : jobs) {
    job.summary.acat_p = job.acat_acc.p();
    if (job.best.p <= pthr && job.best.p <= 1.0)
      write_pair_line(out.top, job.best, Model::Lm, scope);
    if (opt.perm > 0) stage2_perm_topk(opt, Model::Lm, job);
    else {
      job.summary.p_emp = std::numeric_limits<double>::quiet_NaN();
      job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
    }
  }
}

// Explicit instantiations
template void scan_lm_snp_outer<PlinkBed>(const Options&, PlinkBed&, const MissPolicy&, double,
                                          const std::string&, ScopeOut&, double, std::vector<GeneLmJob>&);
template void scan_lm_snp_outer<VcfSession>(const Options&, VcfSession&, const MissPolicy&, double,
                                            const std::string&, ScopeOut&, double, std::vector<GeneLmJob>&);

} // namespace eqtl
