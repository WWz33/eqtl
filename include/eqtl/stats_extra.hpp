/* eqtl — ACAT + empirical p */
#pragma once
#include <vector>
#include <cmath>

namespace eqtl {

double acat(const std::vector<double>& pvals);

// p_emp = (1 + count(T_perm >= T_obs)) / (B+1)
double p_emp_count(double T_obs, const std::vector<double>& T_perm);

// Beta approx of min-p / max-r2 style: fit beta to perm p-values (or 1-r2)
// returns p_beta; shapes optional.
void beta_approx_p(const std::vector<double>& perm_min_p, double obs_min_p,
                   double& p_beta, double& shape1, double& shape2);

} // namespace eqtl
