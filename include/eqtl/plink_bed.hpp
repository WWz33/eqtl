/* eqtl — PLINK bed/bim/fam (GCTA/GEMMA-style) */
#pragma once
#include "eqtl/vcf_session.hpp" // SnpRec, MissPolicy
#include "eqtl/options.hpp"
#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include <fstream>

namespace eqtl {

struct BimSite {
  std::string chrom;
  std::string id;
  double cm = 0;
  int64_t pos = 0;
  std::string a1; // effect allele; dosage = count of A1 (PLINK/GCTA/GEMMA)
  std::string a2;
};

// SNP-major PLINK1 bed. Dosage = number of A1 alleles (0/1/2), same as GEMMA ReadFile_bed.
class PlinkBed {
public:
  PlinkBed() = default;
  ~PlinkBed();
  PlinkBed(const PlinkBed&) = delete;
  PlinkBed& operator=(const PlinkBed&) = delete;

  void open(const std::string& prefix); // expects prefix.bed/.bim/.fam
  const std::string& prefix() const { return prefix_; }
  const std::vector<std::string>& samples() const { return samples_; } // IID order = fam
  const std::vector<BimSite>& sites() const { return sites_; }
  size_t n_samples_file() const { return n_file_; }
  size_t n_snps() const { return sites_.size(); }

  // Analysis sample order (subset of fam IIDs); columns into bed row
  void set_sample_order(const std::vector<std::string>& sample_ids);

  // maf_min: keep if maf_min <= MAF <= 1-maf_min among non-missing analysis samples; 0 = off
  void for_each_snp(const MissPolicy& miss, double maf_min,
                    const std::function<bool(const SnpRec&)>& fn);
  void for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                           const MissPolicy& miss, double maf_min,
                           const std::function<bool(const SnpRec&)>& fn);
  std::vector<SnpRec> load_all(const MissPolicy& miss, double maf_min, int max_snps = -1);
  std::vector<SnpRec> load_region(const std::string& chrom, int64_t start, int64_t end,
                                  const MissPolicy& miss, double maf_min);

private:
  std::string prefix_;
  std::ifstream bed_;
  std::vector<std::string> samples_; // fam IID
  std::vector<BimSite> sites_;
  size_t n_file_ = 0;   // fam lines
  size_t bytes_per_snp_ = 0;
  bool snp_major_ = true;
  std::vector<int> sample_col_; // fam index per analysis sample
  std::vector<uint8_t> rowbuf_;

  void read_fam(const std::string& path);
  void read_bim(const std::string& path);
  void open_bed(const std::string& path);
  bool decode_snp(size_t snp_idx, const MissPolicy& miss, double maf_min, SnpRec& out);
};

// Shared MAF gate (GEMMA: maf < thr || maf > 1-thr → drop; thr<=0 → off)
inline bool pass_maf(double maf, double maf_min) {
  if (maf_min <= 0.0) return true;
  return !(maf < maf_min || maf > (1.0 - maf_min));
}

} // namespace eqtl
