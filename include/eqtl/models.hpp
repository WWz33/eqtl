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
  double tss_dist = 0; // nan if gw
  bool has_tss_dist = false;
  // model extras
  double phi = 0;
  bool glm_converged = true;
  bool glmm_converged = true;
  bool has_phi = false;
};

struct GenePrepLm {
  Eigen::VectorXd y; // residualized / working y after X projection setup
  Eigen::MatrixXd X;
  Eigen::MatrixXd XtX_inv;
  int n = 0;
  int p = 0;
  double yty = 0;
};

// Fit null projection once per gene for lm.
GenePrepLm prep_lm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X);
// Test SNP g (length n).
AssocHit test_lm(const GenePrepLm& prep, const Eigen::VectorXd& g);

struct LmmBasis {
  Eigen::MatrixXd Q;
  Eigen::VectorXd lambda;
};

// Eigen-decompose GRM once per analysis (not per gene).
LmmBasis make_lmm_basis(const Eigen::MatrixXd& K);

struct GenePrepLmm {
  Eigen::VectorXd y_til; // Q' y
  Eigen::MatrixXd X_til;
  Eigen::VectorXd lambda;
  Eigen::MatrixXd Q;
  double delta = 1; // used when fast
  bool fast = false;
  int n = 0;
  int p = 0;
};

GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const LmmBasis& basis,
                     bool fast);
// Convenience: decomposes K each call (prefer make_lmm_basis + prep_lmm).
GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& K,
                     bool fast);
AssocHit test_lmm(const GenePrepLmm& prep, const Eigen::VectorXd& g);

struct GenePrepGlm {
  Eigen::VectorXd y; // counts
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

bool looks_like_counts(const Eigen::VectorXd& y, int check_n = 20);

} // namespace eqtl
