#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <atomic>
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
  // Golden section with the carried-over interior point reused: after the
  // initial bracket each step evaluates only the point it does not already
  // have, which halves the REML count for the same 40 shrink steps (80 -> 42
  // evaluations). Reusing a value instead of recomputing it does not reproduce
  // the old rounding exactly — the converged delta can differ in the last
  // bits — so output comparisons for this path need a tolerance, not `==`.
  double m1 = hi - (hi - lo) / phi;
  double m2 = lo + (hi - lo) / phi;
  double f1 = reml_negll(m1, y_til, X_til, lambda, df);
  double f2 = reml_negll(m2, y_til, X_til, lambda, df);
  for (int it = 0; it < 40; ++it) {
    if (f1 < f2) {
      hi = m2;
      m2 = m1;
      f2 = f1;
      m1 = hi - (hi - lo) / phi;
      f1 = reml_negll(m1, y_til, X_til, lambda, df);
    } else {
      lo = m1;
      m1 = m2;
      f1 = f2;
      m2 = lo + (hi - lo) / phi;
      f2 = reml_negll(m2, y_til, X_til, lambda, df);
    }
  }
  double d = 0.5 * (lo + hi);
  double cur = reml_negll(d, y_til, X_til, lambda, df);
  // Newton polish (GEMMA-style) with numeric derivatives; fallback to golden answer.
  for (int it = 0; it < 8; ++it) {
    const double h = std::max(d * 1e-4, 1e-9);
    const double fm = reml_negll(d - h, y_til, X_til, lambda, df);
    const double fp = reml_negll(d + h, y_til, X_til, lambda, df);
    const double g1 = (fp - fm) / (2.0 * h);
    const double g2 = (fp - 2.0 * cur + fm) / (h * h);
    if (!(g2 > 1e-12)) break;
    const double step = g1 / g2;
    const double dn = d - step;
    if (!(dn > 1e-8 && dn < 1e6)) break;
    const double vn = reml_negll(dn, y_til, X_til, lambda, df);
    if (!(vn < cur)) break;
    d = dn;
    cur = vn;
    if (std::abs(step) < d * 1e-8) break;
  }
  return d;
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
  // null weighted RSS (X only) for partial R² and bordered-Schur caches.
  {
    Eigen::MatrixXd A00 = p.X_til.transpose() * p.dinv.asDiagonal() * p.X_til;
    const Eigen::VectorXd Dy_til = p.dinv.asDiagonal() * p.y_til;
    Eigen::VectorXd XtDy0 = p.X_til.transpose() * Dy_til;
    p.y_dy = p.y_til.dot(Dy_til);
    p.ldlt_a00 = Eigen::LDLT<Eigen::MatrixXd>(A00);
    // LDLT does not flag singular-but-factorizable input; verify the solve.
    Eigen::VectorXd chi0_check = p.ldlt_a00.solve(XtDy0);
    if (p.ldlt_a00.info() == Eigen::Success &&
        (A00 * chi0_check - XtDy0).cwiseAbs().maxCoeff() <=
            1e-8 * std::max(1.0, XtDy0.cwiseAbs().maxCoeff())) {
      p.chi0 = chi0_check;
      p.rss_null = p.y_dy - XtDy0.dot(p.chi0);
      if (p.rss_null < 0) p.rss_null = 0;
      p.has_a00 = true;
    } else {
      static std::atomic<int> warned{0};
      if (!warned.exchange(1))
        warn("lmm: rank-deficient covariates for at least one gene; p=NaN for its SNPs");
      p.ok = false;
      p.rss_null = p.y_dy;
      p.has_a00 = false;
    }
  }
  return p;
}

GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& K,
                     bool fast) {
  return prep_lmm(y, X, make_lmm_basis(K), fast);
}

namespace {
// Per-thread scratch for the bordered-Schur test. test_lmm is called from
// OpenMP workers (cis gene loop, both permutation loops) as well as from
// serial code, so the granularity has to be per thread. The vectors are
// written and consumed within one call — do not re-enter test_lmm between the
// first write and the last read.
struct LmmTestScratch {
  Eigen::VectorXd g_til, Dg, a, u;
};
LmmTestScratch& lmm_test_scratch() {
  static thread_local LmmTestScratch s;
  return s;
}
}  // namespace

AssocHit test_lmm(const GenePrepLmm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  if (!prep.ok) {
    h.p = std::numeric_limits<double>::quiet_NaN();
    return h;
  }
  const int df = prep.n - prep.p - 1;
  if (df <= 0) {
    static std::atomic<int> warned_df{0};
    if (!warned_df.exchange(1))
      warn("lmm: n - p - 1 <= 0 for at least one gene (too few samples or too many covariates); "
           "p=1 for its SNPs");
  }

  if (prep.has_a00) {
    // Bordered information-matrix Schur path. The gene-constant blocks —
    // A00 = X_til^T D X_til (factored once into ldlt_a00), chi0 = A00^{-1}
    // X_til^T D y_til, rss_null, and y_dy — are cached at prep. Per SNP we
    // only pay: Q^T g (np), D·g_til (n), a = X_til^T (D g_til) (np), one
    // p-dimensional triangular solve reusing ldlt_a00 (p²), and a few scalar
    // dot products. The slow path's (p+1)×(p+1) product and repeat LDLT are
    // removed; any SNP-related near-singularity is caught via the Schur S.
    LmmTestScratch& s = lmm_test_scratch();
    s.g_til.noalias() = prep.Q.transpose() * g;
    s.Dg.noalias() = prep.dinv.cwiseProduct(s.g_til);
    const double gg = s.Dg.dot(s.g_til);
    if (gg < 1e-12) { h.p = 1.0; return h; }
    s.a.noalias() = prep.X_til.transpose() * s.Dg;      // X_til^T D g_til
    const double yg = s.Dg.dot(prep.y_til);             // y_til^T D g_til
    s.u = prep.ldlt_a00.solve(s.a);                     // A00^{-1} a
    const double aTu = s.a.dot(s.u);
    const double S = gg - aTu;                          // Schur complement
    if (S <= 1e-15) { h.p = 1.0; return h; }
    const double bg = (yg - s.a.dot(prep.chi0)) / S;    // Wald slope
    // RSS_full = rss_null - S · bg² (bordered elimination of the SNP row);
    // rss_null = y_dy − XtDy0·chi0 is precomputed, β_NS = chi0 − bg·u.
    double q = prep.rss_null - S * bg * bg;
    if (q < 0) q = 0;
    const double sigma2 = (df > 0) ? (q / df) : 1.0;
    h.beta = bg;
    h.se = std::sqrt(std::max(sigma2 / S, 0.0));
    h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
    h.p = p_from_t(h.stat, df);
    h.r2 = (prep.rss_null > 1e-15) ? std::max(0.0, 1.0 - q / prep.rss_null) : 0.0;
    return h;
  }

  // Slow path: rebuild full (p+1)×(p+1) information matrix (kept for the
  // rare case where LDLT(A00) failed at prep time or no cache was built).
  const Eigen::VectorXd g_til = prep.Q.transpose() * g;
  Eigen::MatrixXd Xg(prep.n, prep.p + 1);
  Xg.leftCols(prep.p) = prep.X_til;
  Xg.col(prep.p) = g_til;
  const Eigen::VectorXd& dinv = prep.dinv;
  const double g_wss = g_til.dot(dinv.asDiagonal() * g_til);
  if (g_wss < 1e-12) { h.p = 1.0; return h; }
  Eigen::MatrixXd XtDX = Xg.transpose() * dinv.asDiagonal() * Xg;
  Eigen::VectorXd XtDy = Xg.transpose() * (dinv.asDiagonal() * prep.y_til);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtDX);
  if (ldlt.info() != Eigen::Success) { h.p = 1.0; return h; }
  const Eigen::VectorXd beta = ldlt.solve(XtDy);
  h.beta = beta(prep.p);
  double q = prep.y_til.dot(dinv.asDiagonal() * prep.y_til) - XtDy.dot(beta);
  if (q < 0) q = 0;
  const double sigma2 = (df > 0) ? (q / df) : 1.0;
  Eigen::VectorXd e = Eigen::VectorXd::Zero(prep.p + 1);
  e(prep.p) = 1.0;
  const Eigen::VectorXd cov_col = ldlt.solve(e);
  h.se = std::sqrt(std::max(sigma2 * cov_col(prep.p), 0.0));
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = p_from_t(h.stat, df);
  h.r2 = (prep.rss_null > 1e-15) ? std::max(0.0, 1.0 - q / prep.rss_null) : 0.0;
  return h;
}

} // namespace eqtl
