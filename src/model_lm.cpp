#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>

namespace eqtl {

bool looks_like_counts(const Eigen::VectorXd& y, int check_n) {
  const int n = std::min(check_n, static_cast<int>(y.size()));
  for (int i = 0; i < n; ++i) {
    if (y(i) < 0) return false;
    if (std::fabs(y(i) - std::floor(y(i))) > 1e-8) return false;
  }
  return true;
}

GenePrepLm prep_lm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X) {
  GenePrepLm p;
  p.n = static_cast<int>(y.size());
  p.p = static_cast<int>(X.cols());
  p.X = X;
  p.y = y;
  Eigen::MatrixXd XtX = X.transpose() * X;
  p.XtX_inv = XtX.ldlt().solve(Eigen::MatrixXd::Identity(p.p, p.p));
  // residualize y onto X complement for storage of null residual norm
  Eigen::VectorXd beta = p.XtX_inv * (X.transpose() * y);
  Eigen::VectorXd r = y - X * beta;
  p.yty = r.squaredNorm();
  return p;
}

AssocHit test_lm(const GenePrepLm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  // Project g and y off X: g_s = g - X (X'X)^{-1} X'g
  Eigen::VectorXd Xt_g = prep.X.transpose() * g;
  Eigen::VectorXd g_s = g - prep.X * (prep.XtX_inv * Xt_g);
  Eigen::VectorXd Xt_y = prep.X.transpose() * prep.y;
  Eigen::VectorXd y_s = prep.y - prep.X * (prep.XtX_inv * Xt_y);

  const double gtg = g_s.squaredNorm();
  if (gtg < 1e-12) {
    h.p = 1.0;
    return h;
  }
  const double gty = g_s.dot(y_s);
  h.beta = gty / gtg;
  const double df = prep.n - prep.p - 1;
  if (df <= 0) {
    h.p = 1.0;
    return h;
  }
  const double rss = y_s.squaredNorm() - h.beta * gty;
  const double s2 = std::max(rss / df, 0.0);
  h.se = std::sqrt(s2 / gtg);
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  const double tss = y_s.squaredNorm();
  h.r2 = (tss > 0) ? (1.0 - rss / tss) : 0.0;
  return h;
}

} // namespace eqtl
