/* eqtl — association result + model APIs */
#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>

namespace eqtl {

struct AssocHit {
  std::string gene;
  std::string snp;
  std::string chrom;
  int64_t pos = 0;
  std::string ref;
  std::string alt;
  double maf = 0;
  double beta = 0;
  double se = 0;
  double stat = 0;
  double p = 1;
  double r2 = 0;
  int n = 0;
  double tss_dist = 0;
  bool has_tss_dist = false;
  double phi = 0;
  bool glm_converged = true;
  bool glmm_converged = true;
  bool has_phi = false;
};

struct GenePrepLm {
  Eigen::VectorXd y_s; // y residualized on X (Frisch–Waugh)
  Eigen::MatrixXd X;
  Eigen::MatrixXd XtX_inv;
  int n = 0;
  int p = 0;
  double yty = 0; // ||y_s||^2
};

GenePrepLm prep_lm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X);
AssocHit test_lm(const GenePrepLm& prep, const Eigen::VectorXd& g);

struct LmmBasis {
  Eigen::MatrixXd Q;
  Eigen::VectorXd lambda;
};

LmmBasis make_lmm_basis(const Eigen::MatrixXd& K);
void sparsify_grm(Eigen::MatrixXd& K, double abs_thr = 1e-4);

struct GenePrepLmm {
  Eigen::VectorXd y_til;
  Eigen::MatrixXd X_til;
  Eigen::VectorXd lambda;
  Eigen::MatrixXd Q;
  Eigen::VectorXd dinv; // 1/(delta*lambda+1), null-fixed
  double delta = 1;
  double rss_null = 0; // weighted null RSS (covariates only) → partial R² denom
  int n = 0;
  int p = 0;
};

// Null REML for delta on X only; SNP tests use fixed delta + Wald.
GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const LmmBasis& basis,
                     bool fast = false);
GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& K,
                     bool fast = false);
AssocHit test_lmm(const GenePrepLmm& prep, const Eigen::VectorXd& g);

struct GenePrepGlm {
  Eigen::VectorXd y;
  Eigen::MatrixXd X;
  Eigen::VectorXd offset;
  double phi = 1;
  bool fast = false;
  bool converged = true;
  int n = 0;
};

GenePrepGlm prep_glm_nb(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, bool fast);
AssocHit test_glm_nb(const GenePrepGlm& prep, const Eigen::VectorXd& g);

struct GenePrepGlmm {
  Eigen::VectorXd y;
  Eigen::MatrixXd X;
  Eigen::MatrixXd K;
  Eigen::VectorXd offset;
  double sigma2 = 1;
  Eigen::VectorXd u;
  Eigen::VectorXd mu;
  bool fast = false;
  bool converged = true;
  int n = 0;
};

GenePrepGlmm prep_glmm_pois(const Eigen::VectorXd& y, const Eigen::MatrixXd& X,
                            const Eigen::MatrixXd& K, bool fast);
AssocHit test_glmm_pois(const GenePrepGlmm& prep, const Eigen::VectorXd& g);

// Non-negative integers on finite entries (skips non-finite). Checks up to check_n finite values.
bool looks_like_counts(const Eigen::VectorXd& y, int check_n = 20);

} // namespace eqtl
