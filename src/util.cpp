#include "eqtl/util.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <cctype>
#include <limits>
#include <cstdlib>
#include <stdexcept>

namespace eqtl {

void die(const std::string& msg) {
  std::cerr << "[E] " << msg << "\n";
  // throw so OMP parallel regions don't exit() (UB); main catches
  throw std::runtime_error(msg);
}
void warn(const std::string& msg) {
  std::cerr << "[W] " << msg << "\n";
}
void info(const std::string& msg) {
  std::cerr << "[I] " << msg << "\n";
}

std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

std::vector<std::string> split_ws(const std::string& s) {
  std::vector<std::string> out;
  std::istringstream iss(s);
  std::string t;
  while (iss >> t) out.push_back(t);
  return out;
}

std::vector<std::string> split_char(const std::string& s, char c) {
  std::vector<std::string> out;
  std::string cur;
  for (char ch : s) {
    if (ch == c) {
      out.push_back(cur);
      cur.clear();
    } else cur.push_back(ch);
  }
  out.push_back(cur);
  return out;
}

std::vector<std::string> intersect_order(
    const std::vector<std::string>& prefer_order,
    const std::vector<std::string>& other) {
  auto m = index_map(other);
  std::vector<std::string> out;
  out.reserve(prefer_order.size());
  for (const auto& id : prefer_order) {
    if (m.count(id)) out.push_back(id);
  }
  return out;
}


std::string chrom_key(const std::string& chrom) {
  std::string s = chrom;
  // Strip leading chr/Chr/CHR
  if (s.size() >= 3) {
    char c0 = static_cast<char>(std::tolower(static_cast<unsigned char>(s[0])));
    char c1 = static_cast<char>(std::tolower(static_cast<unsigned char>(s[1])));
    char c2 = static_cast<char>(std::tolower(static_cast<unsigned char>(s[2])));
    if (c0 == 'c' && c1 == 'h' && c2 == 'r') s = s.substr(3);
  }
  // Strip leading zeros: 01 -> 1, 020 -> 20
  size_t start = 0;
  while (start + 1 < s.size() && s[start] == '0') ++start;
  if (start > 0) s = s.substr(start);
  return s;
}

bool chrom_equal(const std::string& a, const std::string& b) {
  if (a == b) return true;
  return chrom_key(a) == chrom_key(b);
}

void validate_chrom_names(const std::vector<std::string>& geno_chroms,
                          const std::vector<std::string>& gff_chroms) {
  std::unordered_set<std::string> geno_keys;
  for (const auto& c : geno_chroms) geno_keys.insert(chrom_key(c));
  std::unordered_set<std::string> gff_keys;
  for (const auto& c : gff_chroms) gff_keys.insert(chrom_key(c));
  int n_overlap = 0;
  for (const auto& k : geno_keys)
    if (gff_keys.count(k)) ++n_overlap;
  if (n_overlap == 0) {
    std::string geno_ex, gff_ex;
    for (size_t i = 0; i < geno_chroms.size() && i < 5; ++i) {
      if (i) geno_ex += ", ";
      geno_ex += geno_chroms[i];
    }
    for (size_t i = 0; i < gff_chroms.size() && i < 5; ++i) {
      if (i) gff_ex += ", ";
      gff_ex += gff_chroms[i];
    }
    die("chromosome name mismatch: genotype has [" + geno_ex +
        (geno_chroms.size() > 5 ? ", ..." : "") +
        "] but GFF has [" + gff_ex +
        (gff_chroms.size() > 5 ? ", ..." : "") +
        "]. No overlap after chr-prefix and leading-zero normalization.");
  }
}

double pnorm_two_sided(double z) {
  if (!std::isfinite(z)) return 1.0;
  // P(|Z|>z) = erfc(|z|/sqrt(2)); std::erfc is accurate in the tails
  double p = std::erfc(std::fabs(z) / std::sqrt(2.0));
  if (p < 0) p = 0;
  if (p > 1) p = 1;
  return p;
}

// Standard-normal quantile Φ⁻¹(p) via Acklam 2003 rational approximation
// refined with one Halley iteration using the high-accuracy std::erfc-based
// normal CDF. After polish, max relative error is at machine-precision
// (~1e-15) over p ∈ (0,1), matching scipy.stats.norm.ppf to the last ULP for
// typical phenotype values (verified bit-identical downstream against an
// external scipy-based normalization on the smoke pheno). p outside (0,1)
// clamps to the supported range; NaN -> NaN preserved.
double qnorm_inv(double p) {
  if (!std::isfinite(p)) return std::numeric_limits<double>::quiet_NaN();
  if (p <= 0.0) p = 1e-15;
  else if (p >= 1.0) p = 1.0 - 1e-15;
  const double p_target = p;  // CDF(v) should equal this post-Halley
  // Lower-tail symmetry: compute Acklam in the (0, 0.5] half and negate.
  bool neg = (p > 0.5);
  if (neg) p = 1.0 - p;
  const double a[6] = {-39.6968302866538, 220.9460984245205, -275.9285104469687,
                       138.3577518678900, -30.6647980661472, 2.506628277459239};
  const double b[5] = {-54.4760987982241, 161.5858368580409, -155.6989798598867,
                       66.8013118877197,  -13.28068155278571};
  const double c[6] = {-0.007784894002430293, -0.3223964580411364, -2.4007582771618388,
                       -2.5497325393437340, 4.3746641414649678, 2.9381639826987831};
  const double d[4] = {0.0077846957090414622, 0.32246712907003983,
                       2.4451341371429960, 3.7544086619074161};
  const double pl = 0.02425;
  double q, r, v;
  if (p < pl) {
    q = std::sqrt(-2.0 * std::log(p));
    v = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
        ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1.0);
  } else {
    q = p - 0.5;
    r = q * q;
    v = ((((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q) /
        (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4]) * r + 1.0);
  }
  if (neg) v = -v;
  // Halley polish: f(x) = CDF(x) - p_target, f' = φ(x), f'' = -x φ(x).
  // x ← x − 2 f f' / (2 f'² − f f'') = x − 2 f φ / (2 φ² + f x φ).
  const double cdf_v = 0.5 * std::erfc(-v / std::sqrt(2.0));
  const double phi_v = std::exp(-0.5 * v * v) / std::sqrt(2.0 * M_PI);
  if (phi_v > 1e-300) {
    const double f = cdf_v - p_target;
    const double denom = 2.0 * phi_v * phi_v + f * v * phi_v;
    if (denom != 0.0) v -= 2.0 * f * phi_v / denom;
  }
  return v;
}

// Rank-based Inverse Normal Transform (Blizzard 2010). Average-rank convention
// for ties → matches scipy.stats.rankdata(method='average'), QTLtools --normal,
// fastqtl --rank. Done on the finite subset; non-finite entries preserved
// (so downstream sample filtering sees the same missingness pattern).
void inverse_normal_transform(Eigen::VectorXd& y) {
  const int n = y.size();
  if (n <= 1) return;
  std::vector<std::pair<double, int>> v;
  v.reserve(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (std::isfinite(y[i])) v.emplace_back(y[i], i);
  }
  const int nf = static_cast<int>(v.size());
  if (nf <= 1) return;
  std::sort(v.begin(), v.end(),
            [](const std::pair<double,int>& a, const std::pair<double,int>& b){
              return a.first < b.first;
            });
  std::vector<double> rank(static_cast<size_t>(n), 0.0);
  // ties share the average of their 1-based positions
  int i = 0;
  while (i < nf) {
    int j = i + 1;
    while (j < nf && v[static_cast<size_t>(j)].first == v[static_cast<size_t>(i)].first) ++j;
    const double avg = 0.5 * (i + 1 + j);  // rank_start + rank_end / 2 (ranks 1..nf)
    for (int k = i; k < j; ++k) rank[static_cast<size_t>(v[static_cast<size_t>(k)].second)] = avg;
    i = j;
  }
  for (int idx = 0; idx < n; ++idx) {
    if (!std::isfinite(y[idx])) continue;
    const double q = (rank[static_cast<size_t>(idx)] - 0.5) / static_cast<double>(nf);
    y[idx] = qnorm_inv(q);
  }
}

// regularized incomplete beta Ix(a,b) continued fraction
static double betacf(double a, double b, double x) {
  const int maxit = 200;
  const double eps = 3e-12;
  const double fpmin = 1e-30;
  double qab = a + b;
  double qap = a + 1.0;
  double qam = a - 1.0;
  double c = 1.0;
  double d = 1.0 - qab * x / qap;
  if (std::fabs(d) < fpmin) d = fpmin;
  d = 1.0 / d;
  double h = d;
  for (int m = 1; m <= maxit; ++m) {
    int m2 = 2 * m;
    double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < fpmin) d = fpmin;
    c = 1.0 + aa / c;
    if (std::fabs(c) < fpmin) c = fpmin;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (std::fabs(d) < fpmin) d = fpmin;
    c = 1.0 + aa / c;
    if (std::fabs(c) < fpmin) c = fpmin;
    d = 1.0 / d;
    double del = d * c;
    h *= del;
    if (std::fabs(del - 1.0) < eps) break;
  }
  return h;
}

static double betai(double a, double b, double x) {
  if (x < 0 || x > 1) return std::numeric_limits<double>::quiet_NaN();
  if (x == 0 || x == 1) return x;
  double lbeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
  double front = std::exp(std::log(x) * a + std::log(1 - x) * b - lbeta) / a;
  if (x < (a + 1) / (a + b + 2)) return front * betacf(a, b, x);
  return 1.0 - front * betacf(b, a, 1 - x) * a / b;
}

double beta_cdf(double x, double a, double b) {
  if (a <= 0 || b <= 0) return std::numeric_limits<double>::quiet_NaN();
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  return betai(a, b, x);
}

double p_from_t(double t, double df) {
  if (!std::isfinite(t) || df <= 0) return 1.0;
  double x = df / (df + t * t);
  double p = betai(0.5 * df, 0.5, x);
  return std::min(1.0, std::max(0.0, p));
}


void assert_unique_ids(const std::vector<std::string>& ids, const std::string& what) {
  std::unordered_set<std::string> seen;
  seen.reserve(ids.size() * 2);
  for (const auto& id : ids) {
    if (!seen.insert(id).second)
      die(what + ": duplicate id: " + id);
  }
}

std::vector<int> dedupe_ids(const std::vector<std::string>& ids,
                            const std::function<bool(int, int)>& payload_equal,
                            const std::string& what) {
  std::unordered_map<std::string, int> first;
  first.reserve(ids.size() * 2);
  std::vector<int> keep;
  keep.reserve(ids.size());
  int n_dup_ok = 0;
  for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
    const auto& id = ids[static_cast<size_t>(i)];
    auto it = first.find(id);
    if (it == first.end()) {
      first.emplace(id, i);
      keep.push_back(i);
      continue;
    }
    const int j = it->second;
    if (!payload_equal(j, i)) {
      die(what + ": duplicate id with conflicting values: " + id);
    }
    ++n_dup_ok;
  }
  if (n_dup_ok > 0)
    warn(what + ": " + std::to_string(n_dup_ok) +
         " duplicate id row(s) with identical values (kept first)");
  return keep;
}

} // namespace eqtl
