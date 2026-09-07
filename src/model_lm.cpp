#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <atomic>
#include <cmath>

namespace eqtl {

// check_n ignored: scan all finite entries (correctness over early exit)
bool looks_like_counts(const Eigen::VectorXd& y, int /*check_n*/) {
  int seen = 0;
  for (int i = 0; i < y.size(); ++i) {
    if (!std::isfinite(y(i))) continue;
    ++seen;
    if (y(i) < 0) return false;
    if (std::fabs(y(i) - std::floor(y(i))) > 1e-8) return false;
  }
  return seen > 0;
}

// Eigen LDLT does not reliably report rank deficiency via info(); verify by
// reconstruction instead (scale-invariant, cost is one p×p product).
bool ldlt_inv_ok(const Eigen::MatrixXd& XtX, Eigen::MatrixXd& inv) {
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtX);
  if (ldlt.info() != Eigen::Success) return false;
  const int p = static_cast<int>(XtX.rows());
  inv = ldlt.solve(Eigen::MatrixXd::Identity(p, p));
  return (XtX * inv - Eigen::MatrixXd::Identity(p, p)).cwiseAbs().maxCoeff() < 1e-6;
}

GenePrepLm prep_lm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X) {
  GenePrepLm p;
  p.n = static_cast<int>(y.size());
  p.p = static_cast<int>(X.cols());
  p.X = X;
  const Eigen::MatrixXd XtX = X.transpose() * X;
  if (!ldlt_inv_ok(XtX, p.XtX_inv)) {
    // rank-deficient covariates: no test is valid for this gene; callers turn
    // every hit into p=NaN so the gene lands in region.tsv with n_tested=0.
    static std::atomic<int> warned{0};
    if (!warned.exchange(1))
      warn("lm: rank-deficient covariates for at least one gene; those genes are untestable");
    p.ok = false;
    // keep shapes valid: scan_lm_snp_outer broadcasts y_s into a row-major Ys
    p.y_s = Eigen::VectorXd::Zero(p.n);
    return p;
  }
  if (p.n - p.p - 1 <= 0) {
    static std::atomic<int> warned_df{0};
    if (!warned_df.exchange(1))
      warn("lm: n - p - 1 <= 0 for at least one gene (too few samples or too many covariates); "
           "p=1 for those genes");
    p.ok = false;
    p.y_s = Eigen::VectorXd::Zero(p.n);
    return p;
  }
  const Eigen::VectorXd beta = p.XtX_inv * (X.transpose() * y);
  p.y_s = y - X * beta;
  p.yty = p.y_s.squaredNorm();
  return p;
}

AssocHit test_lm(const GenePrepLm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  if (!prep.ok) {
    h.p = std::numeric_limits<double>::quiet_NaN();
    return h;
  }
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
