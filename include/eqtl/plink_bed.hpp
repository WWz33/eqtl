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
#include <memory>

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
//
// The read-only metadata parsed from .fam/.bim (samples_, sites_,
// chrom_range_, pair_lut_, n_file_) is held via shared_ptr<const ...> so a
// `clone_for_thread()` companion shares the ~tens-of-MB BIM/FAM tables across
// every OpenMP worker without re-parsing them on each thread. Only the per-
// thread mutable decode state (FILE* bed_fp_, setvbuf buffer, block buffer,
// sample_col_ index, snp_reuse_) is private to each instance. The sites_ vector
// is const so decoders cannot mutate shared metadata.
class PlinkBed {
public:
  PlinkBed() = default;
  ~PlinkBed();
  // Construct a per-thread decoder sharing the metadata of `src` (sites,
  // samples, chrom ranges, pair LUT, n_file_) but owning a private bed FILE*
  // and decode buffers. Caller must then set_sample_order(...) before scan.
  // Re-opens the same prefix's .bed. Populates `out`, which must be default-
  // constructed (PlinkBed only supports move-out via destruction; copy/move
  // assignments are deleted to keep the FILE*/buffers strictly per-instance).
  static void clone_for_thread(const PlinkBed& src, PlinkBed& out);

  PlinkBed(const PlinkBed&) = delete;
  PlinkBed& operator=(const PlinkBed&) = delete;

  void open(const std::string& prefix);
  const std::string& prefix() const { return prefix_; }
  std::vector<std::string> samples() const {
    return sites_meta_ ? sites_meta_->samples : std::vector<std::string>{};
  }
  // const access to shared metadata — safe to read across threads, never mutate
  const std::vector<BimSite>& sites() const { return *sites_; }
  size_t n_samples_file() const { return sites_meta_ ? sites_meta_->n_file : 0; }
  size_t n_snps() const { return sites_->size(); }
  std::vector<std::string> chromosomes() const;
  bool has_index() const { return true; } // bed is always random-access

  void set_sample_order(const std::vector<std::string>& sample_ids);

  void for_each_snp(const MissPolicy& miss, double maf_min,
                    const std::function<bool(const SnpRec&)>& fn);
  void for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                           const MissPolicy& miss, double maf_min,
                           const std::function<bool(const SnpRec&)>& fn);
  std::vector<SnpRec> load_all(const MissPolicy& miss, double maf_min, int max_snps = -1);
  std::vector<SnpRec> load_region(const std::string& chrom, int64_t start, int64_t end,
                                  const MissPolicy& miss, double maf_min);

  // thread count nuthook (parallel cis VCF uses per-thread BGZF threads).
  // bed is single-threaded stdio + per-thread FILE*; no thread setting needed.
  void set_threads(int /*n*/) {}

private:
  static constexpr size_t kMinBlockSnps = 4096;
  static constexpr size_t kTargetBlockBytes = 4 << 20; // ~4MB SNP block

  // Shared (across clones), read-only after open: BIM/FAM tables + chrom
  // index + pair LUT. Kept alive by shared_ptr so cloned decoders retain
  // access even after the master PlinkBed is destroyed.
  struct Meta {
    std::vector<BimSite> sites;
    std::vector<std::string> samples;
    std::unordered_map<std::string, std::pair<size_t, size_t>> chrom_range;
    int8_t pair_lut[4]{};
    size_t n_file = 0;
  };
  std::shared_ptr<const Meta> sites_meta_;
  // const-view shortcuts (derived from sites_meta_); set on open/clone
  const std::vector<BimSite>* sites_ = nullptr;
  const std::unordered_map<std::string, std::pair<size_t, size_t>>* chrom_range_ = nullptr;
  const int8_t (*pair_lut_)[4] = nullptr; // → points into sites_meta_

  std::string prefix_;
  FILE* bed_fp_ = nullptr;
  size_t bytes_per_snp_ = 0;
  size_t block_snps_ = kMinBlockSnps;
  std::vector<int> sample_col_;
  std::vector<uint8_t> block_buf_;
  std::vector<char> file_buf_;
  SnpRec snp_reuse_;

  void read_fam(Meta& m, const std::string& path);
  void read_bim(Meta& m, const std::string& path);
  void build_chrom_ranges(Meta& m);
  void open_bed(const std::string& path);
  void init_meta_views();
  bool seek_snp(size_t snp_idx);
  bool decode_row(size_t snp_idx, const uint8_t* row, const MissPolicy& miss, double maf_min,
                  SnpRec& out);
  bool for_each_range(size_t lo, size_t hi, const MissPolicy& miss, double maf_min,
                      const std::function<bool(const SnpRec&)>& fn);
  bool for_each_range_pos(size_t lo, size_t hi, int64_t p0, int64_t p1, const MissPolicy& miss,
                          double maf_min, const std::function<bool(const SnpRec&)>& fn);
};

} // namespace eqtl
