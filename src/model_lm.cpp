#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>

namespace eqtl {

bool looks_like_counts(const Eigen::VectorXd& y, int check_n) {
  int seen = 0;
  for (int i = 0; i < y.size() && seen < check_n; ++i) {
    if (!std::isfinite(y(i))) continue;
    ++seen;
    if (y(i) < 0) return false;
    if (std::fabs(y(i) - std::floor(y(i))) > 1e-8) return false;
  }
  return seen > 0;
}

GenePrepLm prep_lm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X) {
  GenePrepLm p;
  p.n = static_cast<int>(y.size());
  p.p = static_cast<int>(X.cols());
  p.X = X;
  Eigen::MatrixXd XtX = X.transpose() * X;
  p.XtX_inv = XtX.ldlt().solve(Eigen::MatrixXd::Identity(p.p, p.p));
  const Eigen::VectorXd beta = p.XtX_inv * (X.transpose() * y);
  p.y_s = y - X * beta;
  p.yty = p.y_s.squaredNorm();
  return p;
}

AssocHit test_lm(const GenePrepLm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  const Eigen::VectorXd Xt_g = prep.X.transpose() * g;
  const Eigen::VectorXd g_s = g - prep.X * (prep.XtX_inv * Xt_g);

  const double gtg = g_s.squaredNorm();
  if (gtg < 1e-12) {
    h.p = 1.0;
    return h;
  }
  const double gty = g_s.dot(prep.y_s);
  h.beta = gty / gtg;
  const double df = prep.n - prep.p - 1;
  if (df <= 0) {
    h.p = 1.0;
    return h;
  }
  const double rss = prep.yty - h.beta * gty;
  const double s2 = std::max(rss / df, 0.0);
  h.se = std::sqrt(s2 / gtg);
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  h.r2 = (prep.yty > 0) ? (1.0 - rss / prep.yty) : 0.0;
  return h;
}

} // namespace eqtl
