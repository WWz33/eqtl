#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>

namespace eqtl {

// Poisson GLMM with GRM random effect (PQL). SNP is fixed effect.

namespace {

void pql_fit(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& K,
             const Eigen::VectorXd& offset, Eigen::VectorXd& beta, Eigen::VectorXd& u,
             Eigen::VectorXd& mu, double& sigma2, bool update_sigma2, bool& converged) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  if (beta.size() != p) beta = Eigen::VectorXd::Zero(p);
  if (u.size() != n) u = Eigen::VectorXd::Zero(n);
  if (sigma2 < 0) sigma2 = 0.1;
  converged = false;
  double step = 1.0;

  for (int it = 0; it < 40; ++it) {
    const Eigen::VectorXd eta = X * beta + offset + u;
    mu = eta.array().exp().max(1e-8).matrix();
    const Eigen::VectorXd w = mu;
    Eigen::VectorXd ytil(n);
    for (int i = 0; i < n; ++i) {
      ytil(i) = eta(i) - offset(i) + (y(i) - mu(i)) / mu(i);
    }

    Eigen::MatrixXd Sigma = sigma2 * K;
    for (int i = 0; i < n; ++i) {
      Sigma(i, i) += 1.0 / std::max(w(i), 1e-8);
    }
    const Eigen::LDLT<Eigen::MatrixXd> sldlt(Sigma);
    if (sldlt.info() != Eigen::Success) {
      break;
    }
    const Eigen::MatrixXd SiX = sldlt.solve(X);
    const Eigen::MatrixXd XtSiX = X.transpose() * SiX;
    const Eigen::LDLT<Eigen::MatrixXd> xldlt(XtSiX);
    if (xldlt.info() != Eigen::Success) {
      break;
    }
    const Eigen::VectorXd beta_new = xldlt.solve(X.transpose() * sldlt.solve(ytil));
    const Eigen::VectorXd u_new = sigma2 * K * sldlt.solve(ytil - X * beta_new);

    double sigma2_new = sigma2;
    if (update_sigma2) {
      const Eigen::MatrixXd Si = sldlt.solve(Eigen::MatrixXd::Identity(n, n));
      const Eigen::MatrixXd P = Si - SiX * xldlt.solve(SiX.transpose());
      const Eigen::VectorXd Py = P * ytil;
      const Eigen::VectorXd KPy = K * Py;
      const double score = Py.dot(KPy) - (P * K).trace();
      const double ai = KPy.dot(P * KPy);
      if (ai > 1e-12) {
        sigma2_new = sigma2 + step * (score / ai);
      }
      if (sigma2_new < 0) {
        step *= 0.5;
        sigma2_new = std::max(0.0, sigma2 * 0.5);
      }
    }

    const double diff = std::max((beta_new - beta).cwiseAbs().maxCoeff(),
                                 std::fabs(sigma2_new - sigma2));
    beta = beta_new;
    u = u_new;
    sigma2 = std::max(0.0, sigma2_new);
    if (diff < 1e-5) {
      converged = true;
      break;
    }
  }
  mu = (X * beta + offset + u).array().exp().max(1e-8).matrix();
}

}  // namespace

GenePrepGlmm prep_glmm_pois(const Eigen::VectorXd& y, const Eigen::MatrixXd& X,
                            const Eigen::MatrixXd& K, bool fast) {
  if (!looks_like_counts(y)) {
    die("glmm requires non-negative integer counts");
  }
  GenePrepGlmm prep;
  prep.y = y;
  prep.X = X;
  prep.K = K;
  prep.n = static_cast<int>(y.size());
  prep.fast = fast;
  prep.offset = Eigen::VectorXd::Zero(prep.n);
  prep.sigma2 = 0.1;
  Eigen::VectorXd beta;
  pql_fit(y, X, K, prep.offset, beta, prep.u, prep.mu, prep.sigma2, true, prep.converged);
  return prep;
}

AssocHit test_glmm_pois(const GenePrepGlmm& prep, const Eigen::VectorXd& g) {
  AssocHit hit;
  hit.n = prep.n;

  Eigen::MatrixXd Xg(prep.n, prep.X.cols() + 1);
  Xg.leftCols(prep.X.cols()) = prep.X;
  Xg.col(prep.X.cols()) = g;

  Eigen::VectorXd beta;
  Eigen::VectorXd u = prep.u;
  Eigen::VectorXd mu;
  double sigma2 = prep.sigma2;
  bool conv = false;
  // full: update sigma2; fast: freeze sigma2 from null
  pql_fit(prep.y, Xg, prep.K, prep.offset, beta, u, mu, sigma2, !prep.fast, conv);

  hit.glmm_converged = conv;
  hit.beta = beta(prep.X.cols());

  const Eigen::VectorXd w = mu.cwiseMax(1e-8);
  Eigen::MatrixXd Sigma = sigma2 * prep.K;
  for (int i = 0; i < prep.n; ++i) {
    Sigma(i, i) += 1.0 / w(i);
  }
  const Eigen::LDLT<Eigen::MatrixXd> sldlt(Sigma);
  const Eigen::MatrixXd SiX = sldlt.solve(Xg);
  const Eigen::MatrixXd XtSiX = Xg.transpose() * SiX;
  const Eigen::LDLT<Eigen::MatrixXd> xldlt(XtSiX);
  const Eigen::MatrixXd covb =
      xldlt.solve(Eigen::MatrixXd::Identity(Xg.cols(), Xg.cols()));
  hit.se = std::sqrt(std::max(covb(prep.X.cols(), prep.X.cols()), 0.0));
  hit.stat = (hit.se > 0) ? (hit.beta / hit.se) : 0.0;
  hit.p = pnorm_two_sided(hit.stat);
  hit.r2 = 0;
  return hit;
}

}  // namespace eqtl
