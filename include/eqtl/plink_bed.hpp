/* eqtl — PLINK bed/bim/fam (SNP-major; sequential block I/O) */
#pragma once
#include "eqtl/vcf_session.hpp" // SnpRec, MissPolicy
#include "eqtl/options.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <cstdio>
#include <unordered_map>

namespace eqtl {

struct BimSite {
  std::string chrom;
  std::string id;
  double cm = 0;
  int64_t pos = 0;
  std::string a1; // effect allele; dosage = count of A1
  std::string a2;
};

// SNP-major PLINK1 bed. Dosage = A1 allele count (0/1/2).
class PlinkBed {
public:
  PlinkBed() = default;
  ~PlinkBed();
  PlinkBed(const PlinkBed&) = delete;
  PlinkBed& operator=(const PlinkBed&) = delete;

  void open(const std::string& prefix);
  const std::string& prefix() const { return prefix_; }
  const std::vector<std::string>& samples() const { return samples_; }
  const std::vector<BimSite>& sites() const { return sites_; }
  size_t n_samples_file() const { return n_file_; }
  size_t n_snps() const { return sites_.size(); }

  void set_sample_order(const std::vector<std::string>& sample_ids);

  void for_each_snp(const MissPolicy& miss, double maf_min,
                    const std::function<bool(const SnpRec&)>& fn);
  void for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                           const MissPolicy& miss, double maf_min,
                           const std::function<bool(const SnpRec&)>& fn);
  std::vector<SnpRec> load_all(const MissPolicy& miss, double maf_min, int max_snps = -1);
  std::vector<SnpRec> load_region(const std::string& chrom, int64_t start, int64_t end,
                                  const MissPolicy& miss, double maf_min);

private:
  static constexpr size_t kBlockSnps = 256;

  std::string prefix_;
  FILE* bed_fp_ = nullptr;
  std::vector<std::string> samples_;
  std::vector<BimSite> sites_;
  size_t n_file_ = 0;
  size_t bytes_per_snp_ = 0;
  std::vector<int> sample_col_;
  std::vector<uint8_t> block_buf_; // kBlockSnps * bytes_per_snp_
  // chrom_key -> [lo, hi) in sites_ (PLINK bim is contig-contiguous)
  std::unordered_map<std::string, std::pair<size_t, size_t>> chrom_range_;
  int8_t pair_lut_[4]{}; // 00,01,10,11 → dosage or -1

  void read_fam(const std::string& path);
  void read_bim(const std::string& path);
  void open_bed(const std::string& path);
  void build_chrom_ranges();
  void init_lut();
  bool seek_snp(size_t snp_idx);
  // decode one SNP row already in row (bytes_per_snp_ bytes)
  bool decode_row(size_t snp_idx, const uint8_t* row, const MissPolicy& miss, double maf_min,
                  SnpRec& out);
  // sequential walk [lo, hi)
  bool for_each_range(size_t lo, size_t hi, const MissPolicy& miss, double maf_min,
                      const std::function<bool(const SnpRec&)>& fn);
  // optional pos filter inside range (for region)
  bool for_each_range_pos(size_t lo, size_t hi, int64_t p0, int64_t p1, const MissPolicy& miss,
                          double maf_min, const std::function<bool(const SnpRec&)>& fn);
};

inline bool pass_maf(double maf, double maf_min) {
  if (maf_min <= 0.0) return true;
  return !(maf < maf_min || maf > (1.0 - maf_min));
}

} // namespace eqtl
