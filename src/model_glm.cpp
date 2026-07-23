#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>

namespace eqtl {

// NB GLM, log link: Var = mu + phi*mu^2. --fast: fix phi from null.

static void nb_irls(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::VectorXd& offset,
                    double& phi, Eigen::VectorXd& beta, Eigen::VectorXd& mu, bool estimate_phi,
                    bool& converged) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  beta = Eigen::VectorXd::Zero(p);
  mu = y.cwiseMax(0.1);
  converged = false;
  for (int it = 0; it < 50; ++it) {
    Eigen::VectorXd z(n), w(n);
    for (int i = 0; i < n; ++i) {
      const double m = std::max(mu(i), 1e-8);
      const double var = m + phi * m * m;
      w(i) = (m * m) / std::max(var, 1e-12);
      z(i) = std::log(m) - offset(i) + (y(i) - m) / m;
    }
    const Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) break;
    const Eigen::VectorXd beta_new =
        ldlt.solve(X.transpose() * (w.asDiagonal() * (z - offset)));
    const Eigen::VectorXd mu_new = (X * beta_new + offset).array().exp().matrix();
    const double diff = (beta_new - beta).cwiseAbs().maxCoeff();
    beta = beta_new;
    mu = mu_new;
    if (estimate_phi) {
      double num = 0;
      for (int i = 0; i < n; ++i) {
        const double m = std::max(mu(i), 1e-8);
        num += ((y(i) - m) * (y(i) - m) - m) / (m * m);
      }
      phi = std::max(1e-8, num / static_cast<double>(n));
    }
    if (diff < 1e-6) {
      converged = true;
      break;
    }
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
  return p;
}

AssocHit test_glm_nb(const GenePrepGlm& prep, const Eigen::VectorXd& g) {
  AssocHit h;
  h.n = prep.n;
  h.has_phi = true;
  const double gvar = g.array().square().mean() - std::pow(g.mean(), 2);
  if (gvar < 1e-12) {
    h.p = 1.0;
    return h;
  }
  Eigen::MatrixXd Xg(prep.n, prep.X.cols() + 1);
  Xg.leftCols(prep.X.cols()) = prep.X;
  Xg.col(prep.X.cols()) = g;
  double phi = prep.phi;
  Eigen::VectorXd beta, mu;
  bool conv = false;
  nb_irls(prep.y, Xg, prep.offset, phi, beta, mu, !prep.fast, conv);
  h.phi = phi;
  h.glm_converged = conv;
  h.beta = beta(prep.X.cols());
  Eigen::VectorXd w(prep.n);
  for (int i = 0; i < prep.n; ++i) {
    const double m = std::max(mu(i), 1e-8);
    const double var = m + phi * m * m;
    w(i) = (m * m) / std::max(var, 1e-12);
  }
  const Eigen::MatrixXd XtWX = Xg.transpose() * w.asDiagonal() * Xg;
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    h.p = 1.0;
    return h;
  }
  Eigen::VectorXd e = Eigen::VectorXd::Zero(Xg.cols());
  e(prep.X.cols()) = 1.0;
  const Eigen::VectorXd cov_col = ldlt.solve(e);
  h.se = std::sqrt(std::max(cov_col(prep.X.cols()), 0.0));
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = pnorm_two_sided(h.stat);
  h.r2 = 0;
  return h;
}

} // namespace eqtl
