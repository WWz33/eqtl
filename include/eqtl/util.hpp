/* eqtl — small utils */
#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <cmath>
#include <limits>
#include <cctype>

namespace eqtl {

void die(const std::string& msg);
void warn(const std::string& msg);
void info(const std::string& msg);

std::vector<std::string> split_ws(const std::string& s);
std::vector<std::string> split_char(const std::string& s, char c);
std::string trim(const std::string& s);

// Strip leading chr/Chr/CHR for matching only; empty if all stripped.
std::string chrom_key(const std::string& chrom);
bool chrom_equal(const std::string& a, const std::string& b);

// two-sided normal p from |z|
double pnorm_two_sided(double z);
double p_from_t(double t, double df);

// MAF gate: keep if maf_min <= MAF <= 1-maf_min; maf_min<=0 disables
inline bool pass_maf(double maf, double maf_min) {
  if (maf_min <= 0.0) return true;
  return !(maf < maf_min || maf > (1.0 - maf_min));
}

// beta cdf for beta approx (incomplete beta via continued fraction)
double beta_cdf(double x, double a, double b);

std::vector<std::string> intersect_order(
    const std::vector<std::string>& prefer_order,
    const std::vector<std::string>& other);

template <typename T>
std::unordered_map<std::string, int> index_map(const std::vector<T>& ids) {
  std::unordered_map<std::string, int> m;
  m.reserve(ids.size() * 2);
  for (int i = 0; i < (int)ids.size(); ++i) m[ids[i]] = i;
  return m;
}

} // namespace eqtl
