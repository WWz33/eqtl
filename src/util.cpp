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
