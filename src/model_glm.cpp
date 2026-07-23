#include "eqtl/models.hpp"
#include "eqtl/util.hpp"
#include <cmath>

namespace eqtl {

// Negative binomial GLM with log link.
// Var = mu + phi * mu^2 ; weight for IRLS W = mu / (1 + phi*mu)
// No APL (v1). --fast: estimate phi on null, fix for SNP tests (glm only; LMM always null REML).

static void nb_irls(const Eigen::VectorXd& y, const Eigen::MatrixXd& X,
                    const Eigen::VectorXd& offset, double& phi,
                    Eigen::VectorXd& beta, Eigen::VectorXd& mu,
                    bool estimate_phi, bool& converged) {
  const int n = static_cast<int>(y.size());
  const int p = static_cast<int>(X.cols());
  beta = Eigen::VectorXd::Zero(p);
  // init mu
  mu = y.cwiseMax(0.1);
  converged = false;
  for (int it = 0; it < 50; ++it) {
    Eigen::VectorXd eta = (mu.array().log() - offset.array()).matrix();
    // working response
    Eigen::VectorXd z(n), w(n);
    for (int i = 0; i < n; ++i) {
      double m = std::max(mu(i), 1e-8);
      double var = m + phi * m * m;
      w(i) = (m * m) / std::max(var, 1e-12); // for log link: (dmu/deta)^2 / var = m^2 / var
      z(i) = eta(i) + (y(i) - m) / m;
    }
    Eigen::MatrixXd XtWX = X.transpose() * w.asDiagonal() * X;
    Eigen::VectorXd XtWz = X.transpose() * (w.asDiagonal() * (z - offset));
    // z already on eta scale without offset in some formulations; use z as target for Xb+offset
    XtWz = X.transpose() * (w.asDiagonal() * z) - X.transpose() * (w.asDiagonal() * offset);
    // simpler: regress (z) on X with offset absorbed in z' = z - offset? 
    // standard: eta = Xb + offset; z = eta + (y-mu)/mu; W as above; Xb = (X'WX)^{-1} X'W (z - offset)
    Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
    if (ldlt.info() != Eigen::Success) break;
    Eigen::VectorXd beta_new = ldlt.solve(X.transpose() * (w.asDiagonal() * (z - offset)));
    Eigen::VectorXd eta_new = X * beta_new + offset;
    Eigen::VectorXd mu_new = eta_new.array().exp().matrix();
    double diff = (beta_new - beta).cwiseAbs().maxCoeff();
    beta = beta_new;
    mu = mu_new;
    if (estimate_phi) {
      // method of moments-ish: phi = max(0, mean( ((y-mu)^2 - mu) / mu^2 ))
      double num = 0, den = 0;
      for (int i = 0; i < n; ++i) {
        double m = std::max(mu(i), 1e-8);
        num += ((y(i) - m) * (y(i) - m) - m) / (m * m);
        den += 1.0;
      }
      phi = std::max(1e-8, num / std::max(den, 1.0));
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
  p.offset = y.array().sum() > 0
                 ? Eigen::VectorXd::Constant(p.n, std::log(std::max(y.mean(), 0.1)))
                 : Eigen::VectorXd::Zero(p.n);
  // better offset: log library size proxy = log(sum genes) handled outside; use log(mean) constant
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
  Eigen::MatrixXd Xg(prep.n, prep.X.cols() + 1);
  Xg.leftCols(prep.X.cols()) = prep.X;
  Xg.col(prep.X.cols()) = g;
  double phi = prep.phi;
  Eigen::VectorXd beta, mu;
  bool conv = false;
  bool est_phi = !prep.fast;
  nb_irls(prep.y, Xg, prep.offset, phi, beta, mu, est_phi, conv);
  h.phi = phi;
  h.glm_converged = conv;
  h.beta = beta(prep.X.cols());
  // se from final IRLS W
  Eigen::VectorXd w(prep.n);
  for (int i = 0; i < prep.n; ++i) {
    double m = std::max(mu(i), 1e-8);
    double var = m + phi * m * m;
    w(i) = (m * m) / std::max(var, 1e-12);
  }
  Eigen::MatrixXd XtWX = Xg.transpose() * w.asDiagonal() * Xg;
  Eigen::LDLT<Eigen::MatrixXd> ldlt(XtWX);
  if (ldlt.info() != Eigen::Success) {
    h.p = 1.0;
    return h;
  }
  Eigen::MatrixXd covb = ldlt.solve(Eigen::MatrixXd::Identity(Xg.cols(), Xg.cols()));
  h.se = std::sqrt(std::max(covb(prep.X.cols(), prep.X.cols()), 0.0));
  h.stat = (h.se > 0) ? (h.beta / h.se) : 0.0;
  h.p = pnorm_two_sided(h.stat);
  h.r2 = 0;
  return h;
}

} // namespace eqtl
