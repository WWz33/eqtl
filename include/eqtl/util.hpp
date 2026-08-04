/* eqtl — small utils */
#pragma once
#include <string>
#include <functional>
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

// Validate chromosome names match between genotype and annotation.
// Dies with helpful error if no overlap after chr-prefix normalization.
void validate_chrom_names(const std::vector<std::string>& geno_chroms,
                          const std::vector<std::string>& gff_chroms);

// two-sided normal p from |z|
double pnorm_two_sided(double z);
double p_from_t(double t, double df);

// Stable FNV-1a hash for reproducible perm seeds (std::hash<string> is impl-defined).
inline unsigned fnv1a(const std::string& s) {
  unsigned long long h = 14695981039346656037ULL;
  for (unsigned char c : s) {
    h ^= c;
    h *= 1099511628211ULL;
  }
  return static_cast<unsigned>(h & 0xFFFFFFFFu);
}

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

// Die on any duplicate id (genotype columns / gene names).
void assert_unique_ids(const std::vector<std::string>& ids, const std::string& what);

// Table rows: same id + same payload → warn, keep first; conflicting payload → die.
// Returns keep indices.
std::vector<int> dedupe_ids(const std::vector<std::string>& ids,
                            const std::function<bool(int, int)>& payload_equal,
                            const std::string& what);

template <typename T>
std::unordered_map<std::string, int> index_map(const std::vector<T>& ids) {
  std::unordered_map<std::string, int> m;
  m.reserve(ids.size() * 2);
  for (int i = 0; i < (int)ids.size(); ++i) m[ids[i]] = i;
  return m;
}

} // namespace eqtl
