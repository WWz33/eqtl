#include "eqtl/fission.hpp"
#include "eqtl/pheno.hpp"
#include "eqtl/util.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace eqtl {
namespace {

// ── thinning helpers ────────────────────────────────────────────────────────

static bool is_integer_counts(const Eigen::MatrixXd& Y) {
  for (int j = 0; j < Y.cols(); ++j)
    for (int i = 0; i < Y.rows(); ++i) {
      const double v = Y(i, j);
      if (std::isnan(v)) continue;
      if (v < 0.0 || std::abs(v - std::round(v)) > 1e-9) return false;
    }
  return true;
}

// Binomial thinning: Y1[i,d] ~ Binom(Y[i,d], ε), Y2 = Y - Y1
// Y1 ⊥ Y2 | λ (Poisson/NB property)
static std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
thin_binom(const Eigen::MatrixXd& Y, double eps, unsigned seed) {
  const int n = Y.rows(), G = Y.cols();
  Eigen::MatrixXd Y1(n, G), Y2(n, G);
  std::mt19937 rng(seed);
  for (int j = 0; j < G; ++j)
    for (int i = 0; i < n; ++i) {
      const double v = Y(i, j);
      if (std::isnan(v)) { Y1(i, j) = Y2(i, j) = v; continue; }
      const int cnt = static_cast<int>(std::round(v));
      const int y1 = std::binomial_distribution<int>{cnt, eps}(rng);
      Y1(i, j) = static_cast<double>(y1);
      Y2(i, j) = static_cast<double>(cnt - y1);
    }
  return {Y1, Y2};
}

// Gaussian fission: Y1 = ε·Y + Z, Y2 = (1-ε)·Y - Z
// Z[:,d] ~ N(0, σ²_d · ε·(1-ε))  →  Y1 ⊥ Y2 | μ
static std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
thin_gaussian(const Eigen::MatrixXd& Y, double eps, unsigned seed) {
  const int n = Y.rows(), G = Y.cols();
  Eigen::MatrixXd Y1(n, G), Y2(n, G);

  // per-gene sample variance (NaN-aware)
  std::vector<double> gene_var(static_cast<size_t>(G), 0.0);
  for (int j = 0; j < G; ++j) {
    double mean = 0.0; int cnt = 0;
    for (int i = 0; i < n; ++i)
      if (!std::isnan(Y(i, j))) { mean += Y(i, j); ++cnt; }
    if (cnt < 2) continue;
    mean /= cnt;
    double var = 0.0;
    for (int i = 0; i < n; ++i)
      if (!std::isnan(Y(i, j))) var += (Y(i, j) - mean) * (Y(i, j) - mean);
    gene_var[static_cast<size_t>(j)] = var / (cnt - 1);
  }

  std::mt19937 rng(seed);
  for (int j = 0; j < G; ++j) {
    const double sigma = std::sqrt(gene_var[static_cast<size_t>(j)] * eps * (1.0 - eps));
    std::normal_distribution<double> nd(0.0, sigma > 1e-15 ? sigma : 1e-15);
    for (int i = 0; i < n; ++i) {
      const double v = Y(i, j);
      if (std::isnan(v)) { Y1(i, j) = Y2(i, j) = v; continue; }
      const double z = nd(rng);
      Y1(i, j) = eps * v + z;
      Y2(i, j) = (1.0 - eps) * v - z;
    }
  }
  return {Y1, Y2};
}

// Write matrix TSV compatible with eqtl -e and -c formats
static void write_tsv(const std::string& path, const Eigen::MatrixXd& M,
                      const std::vector<std::string>& row_ids,
                      const std::vector<std::string>& col_ids) {
  std::ofstream f(path);
  if (!f) die("cannot write: " + path);
  f << "sample_id";
  for (const auto& c : col_ids) f << '\t' << c;
  f << '\n';
  for (int i = 0; i < M.rows(); ++i) {
    f << row_ids[static_cast<size_t>(i)];
    for (int j = 0; j < M.cols(); ++j) f << '\t' << M(i, j);
    f << '\n';
  }
}

// ── PEER VB-FA  (Stegle et al. 2012, reimplemented float64) ────────────────
// Model:  Y ≈ X @ W.T + noise
//   Y  : n×G   (samples × genes, mean-centered before entry)
//   W  : G×K   loadings
//   X  : n×K   latent factors
//   Alpha : K  ARD precision per factor (prior on W columns)
//   Eps   : G  per-gene noise precision

struct PeerState {
  Eigen::MatrixXd W;        // G×K
  Eigen::MatrixXd X_E1;     // n×K
  Eigen::MatrixXd X_E2S;    // K×K  = n·cov_X + X.T@X
  Eigen::MatrixXd W_E2S;    // K×K  = Σ_d (Σ_d + w_d·w_d.T)
  Eigen::MatrixXd W_Xprec;  // K×K  = Σ_d eps[d]·(Σ_d + w_d·w_d.T)
  Eigen::VectorXd Alpha_E1; // K    E[alpha_k] = a/b_k
  Eigen::VectorXd Eps_E1;   // G    E[eps_d]   = a/b_d
  // stale copies saved after W update, consumed by Eps update
  Eigen::VectorXd  A_last;     // K   Alpha at W-update time
  Eigen::MatrixXd  XE2S_last;  // K×K X_E2S at W-update time
};

// W update. Fills XtY (K×G) for reuse as YtX.T in Eps update.
static void update_W(PeerState& s, const Eigen::MatrixXd& Y,
                     Eigen::MatrixXd& XtY) {
  const int G = Y.cols(), K = static_cast<int>(s.W.cols());
  s.W_E2S   = Eigen::MatrixXd::Zero(K, K);
  s.W_Xprec = Eigen::MatrixXd::Zero(K, K);

  XtY = s.X_E1.transpose() * Y;  // K×G, one DGEMM shared with Eps update

  const Eigen::MatrixXd diag_A = s.Alpha_E1.asDiagonal();
  const Eigen::MatrixXd I_K    = Eigen::MatrixXd::Identity(K, K);

  for (int d = 0; d < G; ++d) {
    // cov_d = (diag(Alpha) + eps_d · X_E2S)^{-1}  via LDLT
    const Eigen::LDLT<Eigen::MatrixXd> ldlt(diag_A + s.Eps_E1[d] * s.X_E2S);
    const Eigen::VectorXd w   = ldlt.solve(XtY.col(d)) * s.Eps_E1[d];
    const Eigen::MatrixXd cov = ldlt.solve(I_K);
    s.W.row(d) = w.transpose();
    const Eigen::MatrixXd E2 = cov + w * w.transpose();
    s.W_E2S   += E2;
    s.W_Xprec += s.Eps_E1[d] * E2;
  }

  s.A_last    = s.Alpha_E1;  // stale save for Eps update
  s.XE2S_last = s.X_E2S;
}

// Alpha update: Gamma posterior  E[alpha_k] = (pa + G/2) / (pb + W_E2S[k,k]/2)
static void update_Alpha(PeerState& s, int G, double pa, double pb) {
  const double a = pa + 0.5 * G;
  s.Alpha_E1 = a / (pb + 0.5 * s.W_E2S.diagonal().array());
}

// X update: shared covariance; rhs via one DGEMM
static void update_X(PeerState& s, const Eigen::MatrixXd& Y) {
  const int n = Y.rows(), K = static_cast<int>(s.W.cols());
  const Eigen::LDLT<Eigen::MatrixXd> ldlt(s.W_Xprec +
                                           Eigen::MatrixXd::Identity(K, K));
  const Eigen::MatrixXd X_cov = ldlt.solve(Eigen::MatrixXd::Identity(K, K));

  // eps_W[d,k] = Eps[d]·W[d,k] → rhs = Y · eps_W  (one DGEMM n×G×K)
  const Eigen::MatrixXd eps_W = s.W.array().colwise() * s.Eps_E1.array();
  s.X_E1  = Y * eps_W * X_cov;
  s.X_E2S = n * X_cov + s.X_E1.transpose() * s.X_E1;
}

// Eps update. YtX = XtY.T passed in — no extra DGEMM needed.
// Uses current Eps_E1 (previous iteration's) with stale A_last / XE2S_last.
// pheno_E2 = pheno_var + Y.*Y (original PEER includes observation variance 0.01)
static void update_Eps(PeerState& s, const Eigen::MatrixXd& Y,
                       const Eigen::MatrixXd& pheno_E2,  // n×G
                       const Eigen::MatrixXd& YtX,  // G×K = XtY.T
                       double pa, double pb) {
  const int n = Y.rows(), G = Y.cols(), K = static_cast<int>(s.W.cols());
  const double a_eps = pa + 0.5 * n;

  // b1[d] = sum_i pheno_E2[i,d] (includes obs variance)
  const Eigen::VectorXd b1 = pheno_E2.colwise().sum().transpose();

  // b2[d] = YtX[d,:] · W[d,:]  (elementwise then rowwise sum)
  const Eigen::VectorXd b2 = (YtX.array() * s.W.array()).rowwise().sum();

  // b3[d] = tr(X_E2S · (Wcov_stale + w·w.T))
  const Eigen::MatrixXd diag_A = s.A_last.asDiagonal();
  const Eigen::MatrixXd I_K    = Eigen::MatrixXd::Identity(K, K);
  Eigen::VectorXd b3(G);
  for (int d = 0; d < G; ++d) {
    const Eigen::LDLT<Eigen::MatrixXd> ldlt(diag_A + s.XE2S_last * s.Eps_E1[d]);
    const Eigen::MatrixXd Wcov = ldlt.solve(I_K);
    const Eigen::VectorXd w    = s.W.row(d).transpose();
    b3[d] = (s.X_E2S.array() * (Wcov + w * w.transpose()).array()).sum();
  }

  const Eigen::ArrayXd b_eps =
      (pb + 0.5 * (b1.array() - 2.0 * b2.array() + b3.array())).max(1e-10);
  s.Eps_E1 = (a_eps / b_eps).matrix();
}

} // namespace

// ── public API ──────────────────────────────────────────────────────────────

Eigen::MatrixXd peer_factors(const Eigen::MatrixXd& Y_raw, int k,
                              int max_iter, double tol, unsigned init_seed) {
  const int n = Y_raw.rows(), G = Y_raw.cols();
  if (k >= n || k >= G)
    throw std::runtime_error("peer_factors: k must be < min(n_samples, n_genes)");

  // center each gene; NaN → 0 on residual (mean imputation)
  Eigen::MatrixXd Y = Y_raw;
  for (int j = 0; j < G; ++j) {
    double mean = 0.0; int cnt = 0;
    for (int i = 0; i < n; ++i)
      if (!std::isnan(Y(i, j))) { mean += Y(i, j); ++cnt; }
    if (cnt > 0) mean /= cnt;
    for (int i = 0; i < n; ++i)
      Y(i, j) = std::isnan(Y(i, j)) ? 0.0 : Y(i, j) - mean;
  }

  // PEER hyperparameters (defaults from Stegle et al.)
  constexpr double Alpha_pa = 0.001, Alpha_pb = 0.1;
  constexpr double Eps_pa   = 0.1,   Eps_pb   = 10.0;

  PeerState s;
  s.W_E2S    = Eigen::MatrixXd::Zero(k, k);
  s.W_Xprec  = Eigen::MatrixXd::Zero(k, k);
  s.Alpha_E1 = Eigen::VectorXd::Ones(k);
  s.Eps_E1   = Eigen::VectorXd::Constant(G, Eps_pa / Eps_pb);
  s.A_last    = s.Alpha_E1;
  s.XE2S_last = Eigen::MatrixXd::Identity(k, k) * n;

  // init X, W ~ N(0,1) — matches original PEER RANDN init
  {
    std::mt19937 rng(init_seed);
    std::normal_distribution<double> nd;
    s.X_E1 = Eigen::MatrixXd::NullaryExpr(
        n, k, [&](Eigen::Index, Eigen::Index) { return nd(rng); });
    s.W = Eigen::MatrixXd::Zero(G, k);
    for (int j = 0; j < k; ++j)
      for (int i = 0; i < G; ++i)
        s.W(i, j) = nd(rng);
  }
  s.X_E2S = n * Eigen::MatrixXd::Identity(k, k) + s.X_E1.transpose() * s.X_E1;

  // seed Alpha from zero W (tight ARD prior; first iterations release factors)
  update_Alpha(s, G, Alpha_pa, Alpha_pb);

  // pheno_E2 = pheno_var + Y.*Y (original PEER, pheno_var=0.01)
  constexpr double pheno_var = 0.01;
  Eigen::MatrixXd pheno_E2 = Y.array().square() + pheno_var;

  // convergence: original PEER checks |delta_res_var| < var_tolerance every iter
  constexpr double var_tolerance = 1e-8;
  double last_res_var = std::numeric_limits<double>::max();
  Eigen::MatrixXd XtY;

  for (int iter = 0; iter < max_iter; ++iter) {
    update_W(s, Y, XtY);
    update_Alpha(s, G, Alpha_pa, Alpha_pb);
    update_X(s, Y);
    update_Eps(s, Y, pheno_E2, XtY.transpose(), Eps_pa, Eps_pb);

    // convergence: residual variance (matches original PEER var_tolerance check)
    const double res_var =
        (Y - s.X_E1 * s.W.transpose()).array().square().mean();
    if (std::abs(last_res_var - res_var) < var_tolerance) break;
    // ponytail: skip full ELBO calc (original computes it but only uses it
    // when tolerance>0; default tolerance=1e-3 fires on bound, but var_tolerance=1e-8
    // almost always fires first in practice). Matching var_tolerance behavior.
    if (std::abs(last_res_var - res_var) < tol) break;
    last_res_var = res_var;
  }

  return s.X_E1;
}

int run_fission(const Options& opt) {
  PhenoData P = load_pheno(opt.pheno);
  const int n = static_cast<int>(P.sample_ids.size());
  const int G = static_cast<int>(P.gene_ids.size());

  if (opt.peer_factors >= n || opt.peer_factors >= G)
    die("--peer-factors must be < min(n_samples, n_genes)  [got " +
        std::to_string(opt.peer_factors) + ", n=" + std::to_string(n) +
        ", G=" + std::to_string(G) + "]");

  const unsigned seed = static_cast<unsigned>(opt.seed >= 0 ? opt.seed : 42);

  Eigen::MatrixXd Y1, Y2;
  if (is_integer_counts(P.Y)) {
    info("fission: integer counts detected → Binomial thinning");
    auto [a, b] = thin_binom(P.Y, opt.fission_epsilon, seed);
    Y1 = std::move(a); Y2 = std::move(b);
  } else {
    info("fission: continuous data detected → Gaussian fission");
    auto [a, b] = thin_gaussian(P.Y, opt.fission_epsilon, seed);
    Y1 = std::move(a); Y2 = std::move(b);
  }

  info("fission: estimating " + std::to_string(opt.peer_factors) +
       " PEER factors from Y1 (n=" + std::to_string(n) +
       " samples, G=" + std::to_string(G) + " genes)");

  const Eigen::MatrixXd factors =
      peer_factors(Y1, opt.peer_factors, opt.fission_max_iter,
                   opt.fission_tol, seed + 1u);

  const std::string y2_path = opt.out + ".Y2.tsv";
  write_tsv(y2_path, Y2, P.sample_ids, P.gene_ids);
  info("fission: wrote " + y2_path);

  const std::string y1_path = opt.out + ".Y1.tsv";
  write_tsv(y1_path, Y1, P.sample_ids, P.gene_ids);
  info("fission: wrote " + y1_path);

  std::vector<std::string> factor_names;
  factor_names.reserve(static_cast<size_t>(opt.peer_factors));
  for (int k = 0; k < opt.peer_factors; ++k)
    factor_names.push_back("factor" + std::to_string(k + 1));
  const std::string f_path = opt.out + ".factors.tsv";
  write_tsv(f_path, factors, P.sample_ids, factor_names);
  info("fission: wrote " + f_path);

  return 0;
}

} // namespace eqtl
