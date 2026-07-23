#include "eqtl/stats_extra.hpp"
#include "eqtl/util.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace eqtl {

namespace {
constexpr long double kPi = 3.14159265358979323846L;
}

double acat(const std::vector<double>& pvals) {
  if (pvals.empty()) return 1.0;
  long double sum = 0.0L;
  const double n = static_cast<double>(pvals.size());
  int used = 0;
  for (double p : pvals) {
    if (!std::isfinite(p) || p <= 0.0) continue;
    if (p >= 1.0) p = 1.0 - 1.0 / n;
    // ACAT: equal weight tan((0.5-p)*pi); p-value = 0.5 - atan(T)/pi
    sum += std::tan((0.5L - static_cast<long double>(p)) * kPi) / static_cast<long double>(n);
    ++used;
  }
  if (used == 0) return 1.0;
  double ac = 0.5 - static_cast<double>(std::atan(sum) / kPi);
  if (ac < 0) ac = 0;
  if (ac > 1) ac = 1;
  return ac;
}

double p_emp_count(double T_obs, const std::vector<double>& T_perm) {
  int ge = 0;
  for (double t : T_perm) {
    if (t >= T_obs) ++ge;
  }
  const int B = static_cast<int>(T_perm.size());
  return (1.0 + ge) / (1.0 + B);
}

void beta_approx_p(const std::vector<double>& perm_min_p, double obs_min_p, double& p_beta,
                   double& shape1, double& shape2) {
  std::vector<double> x;
  x.reserve(perm_min_p.size());
  for (double v : perm_min_p) {
    if (std::isfinite(v) && v > 0.0 && v < 1.0) x.push_back(v);
  }
  if (x.size() < 10) {
    p_beta = obs_min_p;
    shape1 = shape2 = 1.0;
    return;
  }
  const double mean = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
  double var = 0.0;
  for (double v : x) {
    const double d = v - mean;
    var += d * d;
  }
  var /= static_cast<double>(x.size() - 1);
  if (var <= 0 || mean <= 0 || mean >= 1) {
    p_beta = obs_min_p;
    shape1 = shape2 = 1.0;
    return;
  }
  const double t = mean * (1.0 - mean) / var - 1.0;
  if (t <= 0) {
    p_beta = obs_min_p;
    shape1 = shape2 = 1.0;
    return;
  }
  shape1 = std::max(1e-6, mean * t);
  shape2 = std::max(1e-6, (1.0 - mean) * t);
  double o = obs_min_p;
  if (o < 1e-300) o = 1e-300;
  if (o > 1.0 - 1e-16) o = 1.0 - 1e-16;
  p_beta = beta_cdf(o, shape1, shape2);
  if (!std::isfinite(p_beta)) p_beta = obs_min_p;
}

}  // namespace eqtl
