/* eqtl fission — data fission + PEER factor estimation */
#pragma once
#include "eqtl/options.hpp"
#include <Eigen/Dense>

namespace eqtl {

int run_fission(const Options& opt);

// PEER VB-FA on Y (n×G, float64). Centered internally. Returns X (n×k).
Eigen::MatrixXd peer_factors(const Eigen::MatrixXd& Y, int k,
                              int max_iter = 1000, double tol = 1e-3,
                              unsigned init_seed = 1234);

} // namespace eqtl
