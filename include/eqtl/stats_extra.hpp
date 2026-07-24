/* eqtl — ACAT + empirical p */
#pragma once
#include <vector>
#include <cmath>
#include <limits>

namespace eqtl {

// Online ACAT (equal weight): accumulate tan terms, finalize with used count.
struct AcatAcc {
  long double tan_sum = 0.0L;
  int used = 0;
  void add(double p);
  double p() const;
};

double acat(const std::vector<double>& pvals);

// p_emp = (1 + count(T_perm >= T_obs)) / (B+1)
double p_emp_count(double T_obs, const std::vector<double>& T_perm);

void beta_approx_p(const std::vector<double>& perm_min_p, double obs_min_p,
                   double& p_beta, double& shape1, double& shape2);

} // namespace eqtl
