#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>

namespace eqtl {

// NB GLM, log link: Var = mu + phi*mu^2. --fast: fix phi from null.
// Working response: z = log(mu) + (y-mu)/mu  (no offset in z; offset only in eta = Xb+offset)

static void nb_irls(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::VectorXd& offset,
                    double& phi, Eigen::VectorXd& beta, Eigen::VectorXd& mu, bool estimate_phi,
                    bool& converged) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  const int df = std::max(1, n - p);
  beta = Eigen::VectorXd::Zero(p);
  mu = y.cwiseMax(0.1);
  converged = false;
  for (int it = 0; it < 50; ++it) {
    Eigen::VectorXd z(n), w(n);
    for (int i = 0; i < n; ++i) {
      const double m = std::max(mu(i), 1e-8);
      const double var = m + phi * m * m;
      w(i) = (m * m) / std::max(var, 1e-12);
      // working response on eta scale without offset; design uses eta = Xb + offset
      z(i) = std::log(m) + (y(i) - m) / m;
    }
    // target: E[z] ≈ Xb + offset  ⇒  solve Xb ≈ z - offset
    const Eigen::VectorXd rhs = z - offset;
    const Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) break;
    const Eigen::VectorXd beta_new = ldlt.solve(X.transpose() * (w.asDiagonal() * rhs));
    const Eigen::VectorXd mu_new = (X * beta_new + offset).array().exp().matrix();
    const double diff = (beta_new - beta).cwiseAbs().maxCoeff();
    beta = beta_new;
    mu = mu_new;
    if (diff < 1e-6) {
      converged = true;
      break;
    }
  }
  // ponytail: alternate phi after beta IRLS (not joint inside each iter)
  if (estimate_phi) {
    double num = 0;
    for (int i = 0; i < n; ++i) {
      const double m = std::max(mu(i), 1e-8);
      num += ((y(i) - m) * (y(i) - m) - m) / (m * m);
    }
    phi = std::max(1e-8, num / static_cast<double>(df));
  }
}

GenePrepGlm prep_glm_nb(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, bool fast) {
  if (!looks_like_counts(y)) die("glm (NB) requires non-negative integer counts");
  GenePrepGlm p;
  p.y = y;
  p.X = X;
  p.n = static_cast<int>(y.size());
  p.fast = fast;
  p.offset = Eigen::VectorXd::Zero(p.n);
  p.phi = 1.0;
  Eigen::VectorXd beta, mu;
  nb_irls(y, X, p.offset, p.phi, beta, mu, true, p.converged);
  p.mu = mu;
  // Null working weights and inverse Fisher info, reused by the per-SNP score
  // test so each SNP is one matrix-vector instead of a full IRLS refit.
  p.w.resize(p.n);
  for (int i = 0; i < p.n; ++i) {
    const double m = std::max(mu(i), 1e-8);
    p.w(i) = m / (1.0 + p.phi * m);
  }
  const Eigen::MatrixXd XtWX = X.transpose() * p.w.asDiagonal() * X;
  p.XtWX_inv = XtWX.ldlt().solve(Eigen::MatrixXd::Identity(X.cols(), X.cols()));
  return p;
}

AssocHit test_glm_nb(const GenePrepGlm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  h.has_phi = true;
  h.phi = prep.phi;
  h.glm_converged = prep.converged;
  if (!prep.converged) {
    h.p = std::numeric_limits<double>::quiet_NaN();
    return h;
  }
  const double gvar = g.array().square().mean() - std::pow(g.mean(), 2);
  if (gvar < 1e-12) {
    h.p = 1.0;
    return h;
  }
  // Score test at the null (beta_g = 0), phi fixed at the null estimate:
  //   U   = g^T (y - mu) / (1 + phi*mu)
  //   inf = g^T W g - (X^T W g)^T (X^T W X)^-1 (X^T W g),  W = mu/(1+phi*mu)
  Eigen::VectorXd resid(prep.n);
  for (int i = 0; i < prep.n; ++i)
    resid(i) = (prep.y(i) - prep.mu(i)) / (1.0 + prep.phi * prep.mu(i));
  const double U = g.dot(resid);
  const Eigen::VectorXd XtWg = prep.X.transpose() * (prep.w.asDiagonal() * g);
  const double info = g.dot(prep.w.asDiagonal() * g) - XtWg.dot(prep.XtWX_inv * XtWg);
  if (info < 1e-12) {
    h.p = 1.0;
    return h;
  }
  h.beta = U / info; // one-step estimator from the null
  h.se = 1.0 / std::sqrt(info);
  h.stat = U / std::sqrt(info);
  h.p = pnorm_two_sided(h.stat);
  h.r2 = 0;
  return h;
}

} // namespace eqtl
