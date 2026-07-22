#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>
#include <algorithm>

namespace eqtl {

// Spectral LMM (FaST-LMM / GEMMA-class):
// K = Q diag(lambda) Q'
// rotate: y_til = Q' y, X_til = Q' X
// Var(y_til_j) ∝ delta*lambda_j + 1  (delta = sigma_g^2 / sigma_e^2)
// Estimate delta by 1D grid + golden section on REML profile.

static double reml_negll(double delta, const Eigen::VectorXd& y_til,
                         const Eigen::MatrixXd& X_til,
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
  double sigma2 = q / df;
  double logdet_x = std::log(std::max(XtDX.determinant(), 1e-300));
  // 0.5 * (df*log(sigma2) + log|D| + log|X'D^{-1}X|)
  return 0.5 * (df * std::log(sigma2) + logdet_d + logdet_x);
}

static double optimize_delta(const Eigen::VectorXd& y_til, const Eigen::MatrixXd& X_til,
                             const Eigen::VectorXd& lambda) {
  const int n = static_cast<int>(y_til.size());
  const int p = static_cast<int>(X_til.cols());
  const int df = n - p;
  if (df <= 0) return 1.0;
  // coarse grid then golden section
  double best_d = 1.0;
  double best_ll = reml_negll(1.0, y_til, X_til, lambda, df);
  for (double d = 1e-5; d <= 1e5; d *= 2.0) {
    double v = reml_negll(d, y_til, X_til, lambda, df);
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
    double m1 = hi - (hi - lo) / phi;
    double m2 = lo + (hi - lo) / phi;
    double f1 = reml_negll(m1, y_til, X_til, lambda, df);
    double f2 = reml_negll(m2, y_til, X_til, lambda, df);
    if (f1 < f2) hi = m2;
    else lo = m1;
  }
  return 0.5 * (lo + hi);
}

GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X,
                     const Eigen::MatrixXd& K, bool fast) {
  GenePrepLmm p;
  p.n = static_cast<int>(y.size());
  p.p = static_cast<int>(X.cols());
  p.fast = fast;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
  if (es.info() != Eigen::Success) die("GRM eigen decomposition failed");
  p.Q = es.eigenvectors();
  p.lambda = es.eigenvalues().cwiseMax(0.0);
  p.y_til = p.Q.transpose() * y;
  p.X_til = p.Q.transpose() * X;
  if (fast) {
    p.delta = optimize_delta(p.y_til, p.X_til, p.lambda);
  }
  return p;
}

AssocHit test_lmm(const GenePrepLmm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  Eigen::VectorXd g_til = prep.Q.transpose() * g;

  // design: [X_til | g_til]
  Eigen::MatrixXd Xg(prep.n, prep.p + 1);
  Xg.leftCols(prep.p) = prep.X_til;
  Xg.col(prep.p) = g_til;

  double delta = prep.fast ? prep.delta
                           : optimize_delta(prep.y_til, Xg, prep.lambda);

  Eigen::VectorXd dinv(prep.n);
  for (int i = 0; i < prep.n; ++i) {
    double v = delta * prep.lambda(i) + 1.0;
    if (v < 1e-12) v = 1e-12;
    dinv(i) = 1.0 / v;
  }
  Eigen::MatrixXd XtDX = Xg.transpose() * dinv.asDiagonal() * Xg;
  Eigen::VectorXd XtDy = Xg.transpose() * (dinv.asDiagonal() * prep.y_til);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
  if (ldlt.info() != Eigen::Success) {
    h.p = 1.0;
    return h;
  }
  Eigen::VectorXd beta = ldlt.solve(XtDy);
  h.beta = beta(prep.p);
  const int df = prep.n - prep.p - 1;
  double q = prep.y_til.dot(dinv.asDiagonal() * prep.y_til) - XtDy.dot(beta);
  if (q < 0) q = 0;
  double sigma2 = (df > 0) ? (q / df) : 1.0;
  Eigen::MatrixXd covb = sigma2 * ldlt.solve(Eigen::MatrixXd::Identity(prep.p + 1, prep.p + 1));
  h.se = std::sqrt(std::max(covb(prep.p, prep.p), 0.0));
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  // r2 proxy on rotated residual scale
  Eigen::VectorXd fit = Xg * beta;
  Eigen::VectorXd r = prep.y_til - fit;
  double tss = prep.y_til.squaredNorm();
  h.r2 = (tss > 0) ? (1.0 - r.squaredNorm() / tss) : 0.0;
  return h;
}

} // namespace eqtl
