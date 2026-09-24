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
  bool ok = true; // false ⇒ X rank-deficient for this gene (tests return NaN p)
};

// LDLT has no reliable rank-deficiency reporting; verify the computed inverse by
// reconstruction (scale-invariant). Returns false when X'X is not invertible.
bool ldlt_inv_ok(const Eigen::MatrixXd& XtX, Eigen::MatrixXd& inv);

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
  // Per-SNP Wald hot-path caches. The bordered information matrix is
  //   M = [ A00   a ]      a   = X_til^T D g_til   (p)
  //       [ a^T  gg ],     gg  = g_til^T D g_til   (scalar)
  // with constant gene blocks A00 = X_til^T D X_til (p×p), y_dy =
  // y_til^T D y_til (scalar), XtDy0 = X_til^T D y_til (p), and precomputed
  // chi0 = A00^{-1} XtDy0 (p). LDLT of A00 is factored once and reused at
  // every SNP via the bordered Schur complement, eliminating the per-SNP
  // (p+1)×(p+1) product and factorization. has_a00 is false when the LDLT
  // failed at prep time, in which case test_lmm falls back to the slow path.
  bool has_a00 = false;
  bool ok = true; // false ⇒ covariates rank-deficient for this gene (tests return NaN p)
  Eigen::LDLT<Eigen::MatrixXd> ldlt_a00{Eigen::MatrixXd(0,0)};
  Eigen::VectorXd chi0;       // A00^{-1} XtDy0
  double y_dy = 0.0;          // y_til^T D y_til
};

// Re-prep hints for a gene whose y changed but whose model did not (the
// permutation loops): X_til = Q^T X ignores y, and a caller that tests in the
// spectral domain never reads prep.Q, so neither has to be recomputed or
// stored per draw.
// fixed_delta skips the REML search and holds the null at that delta, which
// changes the permutation null slightly: opt-in, see --perm-freeze-delta.
// Contract: x_til, when set, is Q^T X for the same X and basis this prep will
// use. keep_q=false leaves the resulting prep.Q empty, so such a prep may only
// be consumed by tests that work from g_til (test_lmm_gtil).
struct LmmPrepReuse {
  const Eigen::MatrixXd* x_til = nullptr;  // reuse this instead of Q^T X
  bool keep_q = true;                      // false: leave prep.Q empty
  const double* fixed_delta = nullptr;     // skip the REML search, use this
};

// Null REML for delta on X only; SNP tests use fixed delta + Wald.
GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const LmmBasis& basis,
                     bool fast = false, const LmmPrepReuse* reuse = nullptr);
GenePrepLmm prep_lmm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, const Eigen::MatrixXd& K,
                     bool fast = false);
AssocHit test_lmm(const GenePrepLmm& prep, const Eigen::VectorXd& g);

struct GenePrepGlm {
  Eigen::VectorXd y;
  Eigen::MatrixXd X;
  Eigen::VectorXd offset;
  Eigen::VectorXd mu; // null fitted mean (for parametric permutation)
  Eigen::VectorXd w;  // null working weights mu/(1+phi*mu)
  Eigen::MatrixXd XtWX_inv; // inverse null Fisher info for covariates
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
