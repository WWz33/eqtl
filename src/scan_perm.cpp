#include "eqtl/scan_trans.hpp"
#include <random>
#include <atomic>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>

namespace eqtl {

template <typename Job>
void stage2_perm_topk(const Options& opt, Model model, Job& job,
                      const LmmBasis* ext_basis) {
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
  grb.K_ref = job.gr.K_ref ? job.gr.K_ref : &job.gr.K;
  if (ext_basis) {
    grb.basis_ref = ext_basis;
    grb.has_basis = true;
  } else {
    grb.basis_ref = &job.gr.basis;
    grb.has_basis = job.gr.has_basis;
  }
  grb.y = job.gr.y;

  GenePrepLm lm_c;
  GenePrepLmm lmm_c;
  GenePrepGlm glm_c;
  GenePrepGlmm glmm_c;
  prep_null(model, opt.fast, grb, &lm_c, &lmm_c, &glm_c, &glmm_c);

  Eigen::VectorXd y_perm_base = grb.y;
  Eigen::VectorXd Xb0_til;   // LMM only: precomputed null-fitted spectral mean
  if (model == Model::Lm && lm_c.n > 0) {
    y_perm_base = lm_c.y_s;
  } else if (model == Model::Lmm && lmm_c.n > 0) {
    // whitened spectral residuals: Var(w_i) = σ² under null → exchangeable
    const Eigen::VectorXd& dinv = lmm_c.dinv;
    if (lmm_c.has_a00) {
      // prep_lmm has already cached chi0 = A00^{-1} X_til^T D y_til; reuse it
      // instead of re-forming and re-factoring A00 here.
      Xb0_til = lmm_c.X_til * lmm_c.chi0;
      const Eigen::VectorXd r_til = lmm_c.y_til - Xb0_til;
      y_perm_base = r_til.cwiseProduct(dinv.cwiseSqrt());
    } else {
      Eigen::MatrixXd XtDX = lmm_c.X_til.transpose() * dinv.asDiagonal() * lmm_c.X_til;
      Eigen::VectorXd XtDy = lmm_c.X_til.transpose() * (dinv.asDiagonal() * lmm_c.y_til);
      Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
      if (ldlt.info() == Eigen::Success) {
        const Eigen::VectorXd b0 = ldlt.solve(XtDy);
        Xb0_til = lmm_c.X_til * b0;
        const Eigen::VectorXd r_til = lmm_c.y_til - Xb0_til;
        y_perm_base = r_til.cwiseProduct(dinv.cwiseSqrt());
      }
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
                         static_cast<unsigned>(b) * 9973u + fnv1a(job.gene));
      std::vector<int> idx(static_cast<size_t>(y_perm_base.size()));
      std::iota(idx.begin(), idx.end(), 0);
      std::shuffle(idx.begin(), idx.end(), rng_b);
      for (int i = 0; i < y_perm_base.size(); ++i)
        grb2.y(i) = y_perm_base(idx[static_cast<size_t>(i)]);
      if (model == Model::Lmm && lmm_c.n > 0 && lmm_c.Q.size() > 0) {
        // b0 = A00^{-1} X_til^T D y_til (permutation-invariant, cached as
        // lmm_c.chi0; X_til * b0 cached once as Xb0_til before this loop).
        // Re-deriving it inside the loop was pure waste.
        const Eigen::VectorXd& dinv = lmm_c.dinv;
        if (lmm_c.has_a00 && Xb0_til.size() == dinv.size()) {
          // grb2.y holds shuffled whitened residuals w_perm → un-whiten to
          // spectral residuals, add the null-fitted spectral mean, project back
          // to the sample domain so prep_null (on grb2) treats it as a y.
          Eigen::VectorXd r_til_perm(grb2.y.size());
          for (int i = 0; i < grb2.y.size(); ++i)
            r_til_perm(i) = grb2.y(i) / std::sqrt(std::max(dinv(i), 1e-12));
          grb2.y = lmm_c.Q * (Xb0_til + r_til_perm);
        } else {
          Eigen::MatrixXd XtDX = lmm_c.X_til.transpose() * dinv.asDiagonal() * lmm_c.X_til;
          Eigen::VectorXd XtDy = lmm_c.X_til.transpose() * (dinv.asDiagonal() * lmm_c.y_til);
          Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
          if (ldlt.info() == Eigen::Success) {
            const Eigen::VectorXd b0 = ldlt.solve(XtDy);
            Eigen::VectorXd r_til_perm(grb2.y.size());
            for (int i = 0; i < grb2.y.size(); ++i)
              r_til_perm(i) = grb2.y(i) / std::sqrt(std::max(dinv(i), 1e-12));
            grb2.y = lmm_c.Q * (lmm_c.X_til * b0 + r_til_perm);
          }
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

// Explicit instantiations
template void stage2_perm_topk<GeneLmJob>(const Options&, Model, GeneLmJob&, const LmmBasis*);
template void stage2_perm_topk<GeneLmmJob>(const Options&, Model, GeneLmmJob&, const LmmBasis*);

} // namespace eqtl
