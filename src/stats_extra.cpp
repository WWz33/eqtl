#include "eqtl/stats_extra.hpp"
#include "eqtl/util.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

namespace eqtl {

namespace {
constexpr long double kPi = 3.14159265358979323846L;
}

void AcatAcc::add(double p) {
  if (!std::isfinite(p) || p < 0.0) return;
  if (p >= 1.0 - 1e-12) return; // uninformative (degenerate test); would dominate tan sum
  if (p <= 0.0) p = 1e-300; // keep strongest signals
  tan_sum += std::tan((0.5L - static_cast<long double>(p)) * kPi);
  ++used;
}

double AcatAcc::p() const {
  if (used == 0) return 1.0;
  const long double T = tan_sum / static_cast<long double>(used);
  double ac = 0.5 - static_cast<double>(std::atan(T) / kPi);
  if (ac < 0) ac = 0;
  if (ac > 1) ac = 1;
  return ac;
}

double acat(const std::vector<double>& pvals) {
  AcatAcc a;
  for (double p : pvals) a.add(p);
  return a.p();
}

double p_emp_count(double T_obs, const std::vector<double>& T_perm) {
  int ge = 0;
  for (double t : T_perm) {
    if (t >= T_obs) ++ge;
  }
  const int B = static_cast<int>(T_perm.size());
  return (1.0 + ge) / (1.0 + B);
}

namespace {
constexpr double kS1Min = 0.1, kS1Max = 10.0, kS2Min = 5.0, kS2Max = 1e6;

inline bool in_shape_bounds(double a, double b) {
  return a >= kS1Min && a <= kS1Max && b >= kS2Min && b <= kS2Max;
}

struct BetaStats { double s1 = 0, s2 = 0, n = 0; };  // sum log p, sum log(1-p), count

BetaStats beta_sufficient_stats(std::vector<double>& pval) {
  BetaStats st;
  for (double& v : pval) {
    if (v >= 1.0) v = 0.99999999;  // QTLtools convention
    st.s1 += std::log(v);
    st.s2 += std::log(1.0 - v);
    ++st.n;
  }
  return st;
}

double beta_nll(double a, double b, const BetaStats& st) {
  return -((a - 1.0) * st.s1 + (b - 1.0) * st.s2 - st.n * gsl_sf_lnbeta(a, b));
}

double beta_simplex_nll(const gsl_vector* v, void* params) {
  const double a = gsl_vector_get(v, 0), b = gsl_vector_get(v, 1);
  if (!in_shape_bounds(a, b)) return 1e300;
  return beta_nll(a, b, *static_cast<BetaStats*>(params));
}
}  // namespace

bool beta_fit_mm(std::vector<double>& pval, double& shape1, double& shape2) {
  if (pval.size() < 2) return false;
  const double mean = std::accumulate(pval.begin(), pval.end(), 0.0) / pval.size();
  double var = 0.0;
  for (double v : pval) { const double d = v - mean; var += d * d; }
  var /= static_cast<double>(pval.size() - 1);
  if (var <= 0 || mean <= 0 || mean >= 1) return false;
  const double t = mean * (1.0 - mean) / var - 1.0;
  if (t <= 0) return false;
  shape1 = mean * t;
  shape2 = (1.0 - mean) * t;
  return true;
}

// Faithful port of QTLtools learnBetaParameters (cis_learn_beta.cpp):
// GSL Nelder-Mead from MM start, max 1000 iterations, size tolerance 0.01.
bool beta_fit_ml_simplex(std::vector<double>& pval, double& shape1, double& shape2) {
  BetaStats st = beta_sufficient_stats(pval);
  if (shape1 <= 0 || shape2 <= 0) return false;

  gsl_vector* x = gsl_vector_alloc(2);
  gsl_vector_set(x, 0, shape1);
  gsl_vector_set(x, 1, shape2);
  gsl_vector* ss = gsl_vector_alloc(2);
  gsl_vector_set(ss, 0, shape1 / 10.0);
  gsl_vector_set(ss, 1, shape2 / 10.0);
  gsl_multimin_function f;
  f.n = 2;
  f.f = beta_simplex_nll;
  f.params = &st;
  gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 2);
  gsl_multimin_fminimizer_set(s, &f, x, ss);
  size_t iter = 0;
  int status;
  do {
    ++iter;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status) break;
    status = gsl_multimin_test_size(gsl_multimin_fminimizer_size(s), 0.01);
  } while (status == GSL_CONTINUE && iter < 1000);
  shape1 = gsl_vector_get(s->x, 0);
  shape2 = gsl_vector_get(s->x, 1);
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free(s);
  return status == GSL_SUCCESS && std::isfinite(shape1) && std::isfinite(shape2) &&
         in_shape_bounds(shape1, shape2);
}

// Damped Newton was benchmarked against this simplex port (3000 simulated genes):
// identical optimum (max rel diff 2.7e-4), only 1.2x faster (15.4 vs 18.7us) — not worth it.

void beta_approx_p(const std::vector<double>& perm_min_p, double obs_min_p, double& p_beta,
                   double& shape1, double& shape2) {
  std::vector<double> x;
  x.reserve(perm_min_p.size());
  for (double v : perm_min_p) {
    if (std::isfinite(v) && v > 0.0 && v < 1.0) x.push_back(v);
  }
  if (x.size() < 10) {
    p_beta = std::numeric_limits<double>::quiet_NaN();
    shape1 = shape2 = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  if (!beta_fit_mm(x, shape1, shape2)) {
    p_beta = std::numeric_limits<double>::quiet_NaN();
    shape1 = shape2 = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  const double mm1 = shape1, mm2 = shape2;
  if (!beta_fit_ml_simplex(x, shape1, shape2)) {
    shape1 = mm1;  // ML failed: fall back to moment matching (QTLtools behavior)
    shape2 = mm2;
  }
  double o = obs_min_p;
  if (o < 1e-300) o = 1e-300;
  if (o > 1.0 - 1e-16) o = 1.0 - 1e-16;
  p_beta = beta_cdf(o, shape1, shape2);
  if (!std::isfinite(p_beta)) p_beta = std::numeric_limits<double>::quiet_NaN();
}

}  // namespace eqtl
