#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>
#include <algorithm>

namespace eqtl {

// Spectral LMM: K=QΛQ', Var∝δλ+1; null REML on X then fixed-δ Wald. No per-SNP re-REML.

static double reml_negll(double delta, const Eigen::VectorXd& y_til, const Eigen::MatrixXd& X_til,
                         const Eigen::VectorXd& lambda, int df) {
  const int n = static_cast<int>(y_til.size());
  Eigen::VectorXd dinv(n);
  double logdet_d = 0;
  for (int i = 0; i < n; ++i) {
    double v = delta * lambda(i) + 1.0;
    if (v < 1e-12) v = 1e-12;
    dinv(i) = 1.0 / v;
    logdet_d += std::log(v);
  }
  Eigen::MatrixXd XtDX = X_til.transpose() * dinv.asDiagonal() * X_til;
  Eigen::VectorXd XtDy = X_til.transpose() * (dinv.asDiagonal() * y_til);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
  if (ldlt.info() != Eigen::Success) return 1e300;
  Eigen::VectorXd beta = ldlt.solve(XtDy);
  double q = y_til.dot(dinv.asDiagonal() * y_til) - XtDy.dot(beta);
  if (q <= 0) q = 1e-12;
  const double sigma2 = q / df;
  double logdet_x = 0.0;
  const auto& D = ldlt.vectorD();
  for (int i = 0; i < D.size(); ++i) {
    const double di = D(i);
    if (di <= 0) return 1e300;
    logdet_x += std::log(di);
  }
  return 0.5 * (df * std::log(sigma2) + logdet_d + logdet_x);
}

static double optimize_delta(const Eigen::VectorXd& y_til, const Eigen::MatrixXd& X_til,
                             const Eigen::VectorXd& lambda) {
  const int n = static_cast<int>(y_til.size());
  const int p = static_cast<int>(X_til.cols());
  const int df = n - p;
  if (df <= 0) return 1.0;
  double best_d = 1.0;
  double best_ll = reml_negll(1.0, y_til, X_til, lambda, df);
  for (double d = 1e-5; d <= 1e5; d *= 2.0) {
    const double v = reml_negll(d, y_til, X_til, lambda, df);
    if (v < best_ll) {
      best_ll = v;
      best_d = d;
    }
  }
  double lo = best_d / 4.0, hi = best_d * 4.0;
  if (lo < 1e-8) lo = 1e-8;
  if (hi > 1e6) hi = 1e6;
  const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
  for (int it = 0; it < 40; ++it) {
    const double m1 = hi - (hi - lo) / phi;
    const double m2 = lo + (hi - lo) / phi;
    const double f1 = reml_negll(m1, y_til, X_til, lambda, df);
    const double f2 = reml_negll(m2, y_til, X_til, lambda, df);
    if (f1 < f2) hi = m2;
    else lo = m1;
  }
  return 0.5 * (lo + hi);
}

void sparsify_grm(Eigen::MatrixXd& K, double abs_thr) {
  if (abs_thr <= 0.0) return;
  const int n = static_cast<int>(K.rows());
  if (K.cols() != n) return;
  size_t n_zero = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (std::abs(K(i, j)) < abs_thr) {
        K(i, j) = 0.0;
        K(j, i) = 0.0;
        n_zero += 2;
      }
    }
  }
  info("fast: GRM sparse approx thr=" + std::to_string(abs_thr) + " zeroed " +
       std::to_string(n_zero) + " off-diagonal entries");
}

LmmBasis make_lmm_basis(const Eigen::MatrixXd& K) {
  LmmBasis b;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
  if (es.info() != Eigen::Success) die("GRM eigen decomposition failed");
  b.Q = es.eigenvectors();
  b.lambda = es.eigenvalues().cwiseMax(0.0);
  return b;
}

static void fill_dinv(GenePrepLmm& p) {
  p.dinv.resize(p.n);
  for (int i = 0; i < p.n; ++i) {
    double v = p.delta * p.lambda(i) + 1.0;
    if (v < 1e-12) v = 1e-12;
    p.dinv(i) = 1.0 / v;
  }
}

GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const LmmBasis& basis,
                     bool /*fast*/) {
  GenePrepLmm p;
  p.n = static_cast<int>(y.size());
  p.p = static_cast<int>(X.cols());
  p.Q = basis.Q;
  p.lambda = basis.lambda;
  p.y_til = p.Q.transpose() * y;
  p.X_til = p.Q.transpose() * X;
  p.delta = optimize_delta(p.y_til, p.X_til, p.lambda);
  fill_dinv(p);
  // null weighted RSS (X only) for partial R²
  {
    Eigen::MatrixXd XtDX0 = p.X_til.transpose() * p.dinv.asDiagonal() * p.X_til;
    Eigen::VectorXd XtDy0 = p.X_til.transpose() * (p.dinv.asDiagonal() * p.y_til);
    Eigen::LDLT<Eigen::MatrixXd> ldlt0(XtDX0);
    if (ldlt0.info() == Eigen::Success) {
      const Eigen::VectorXd b0 = ldlt0.solve(XtDy0);
      p.rss_null = p.y_til.dot(p.dinv.asDiagonal() * p.y_til) - XtDy0.dot(b0);
      if (p.rss_null < 0) p.rss_null = 0;
    } else {
      p.rss_null = p.y_til.dot(p.dinv.asDiagonal() * p.y_til);
    }
  }
  return p;
}

GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& K,
                     bool fast) {
  return prep_lmm(y, X, make_lmm_basis(K), fast);
}

AssocHit test_lmm(const GenePrepLmm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  const Eigen::VectorXd g_til = prep.Q.transpose() * g;

  Eigen::MatrixXd Xg(prep.n, prep.p + 1);
  Xg.leftCols(prep.p) = prep.X_til;
  Xg.col(prep.p) = g_til;

  const Eigen::VectorXd& dinv = prep.dinv;
  // near-null SNP variance in spectral space → skip unstable fit
  const double g_wss = g_til.dot(dinv.asDiagonal() * g_til);
  if (g_wss < 1e-12) {
    h.p = 1.0;
    return h;
  }
  Eigen::MatrixXd XtDX = Xg.transpose() * dinv.asDiagonal() * Xg;
  Eigen::VectorXd XtDy = Xg.transpose() * (dinv.asDiagonal() * prep.y_til);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
  if (ldlt.info() != Eigen::Success) {
    h.p = 1.0;
    return h;
  }
  const Eigen::VectorXd beta = ldlt.solve(XtDy);
  h.beta = beta(prep.p);
  const int df = prep.n - prep.p - 1;
  double q = prep.y_til.dot(dinv.asDiagonal() * prep.y_til) - XtDy.dot(beta);
  if (q < 0) q = 0;
  const double sigma2 = (df > 0) ? (q / df) : 1.0;
  Eigen::VectorXd e = Eigen::VectorXd::Zero(prep.p + 1);
  e(prep.p) = 1.0;
  const Eigen::VectorXd cov_col = ldlt.solve(e);
  h.se = std::sqrt(std::max(sigma2 * cov_col(prep.p), 0.0));
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  // partial R²: 1 - RSS_full / RSS_null (null = covariates only; matches LM y_s style)
  h.r2 = (prep.rss_null > 1e-15) ? std::max(0.0, 1.0 - q / prep.rss_null) : 0.0;
  return h;
}

} // namespace eqtl
