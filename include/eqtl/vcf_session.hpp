/* eqtl — VCF/BCF session, GT dosage */
#pragma once
#include <string>
#include <vector>
#include <cstdint>
#include <functional>
#include "eqtl/options.hpp"

namespace eqtl {

struct SnpRec {
  std::string chrom;
  int64_t pos = 0;
  std::string id;
  std::string ref;
  std::string alt;
  std::vector<double> dosage; // n samples
  double maf = 0;
};

struct MissPolicy {
  MissHand hand = MissHand::Filter;
  double max_miss = 0.0; // drop if n_miss/n > max_miss
};

class VcfSession {
public:
  VcfSession() = default;
  ~VcfSession();
  VcfSession(const VcfSession&) = delete;
  VcfSession& operator=(const VcfSession&) = delete;

  void open(const std::string& path);
  const std::string& path() const { return path_; }
  const std::vector<std::string>& samples() const { return samples_; }
  const std::vector<std::string>& contigs() const { return contigs_; }
  bool has_index() const { return indexed_; }

  void set_sample_order(const std::vector<std::string>& sample_ids);

  void for_each_snp(const MissPolicy& miss, const std::function<bool(const SnpRec&)>& fn);
  // 1-based inclusive region; uses CSI/TBI when present
  void for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                           const MissPolicy& miss, const std::function<bool(const SnpRec&)>& fn);
  std::vector<SnpRec> load_all(const MissPolicy& miss, int max_snps = -1);
  std::vector<SnpRec> load_region(const std::string& chrom, int64_t start, int64_t end,
                                  const MissPolicy& miss);

private:
  std::string path_;
  void* fp_ = nullptr;
  void* hdr_ = nullptr;
  void* active_hdr_ = nullptr;
  void* tbx_ = nullptr;
  void* idx_ = nullptr;
  bool indexed_ = false;
  bool warned_no_index_ = false;
  std::vector<std::string> samples_;
  std::vector<std::string> contigs_;
  std::vector<int> sample_col_;

  bool parse_record(void* rec, const MissPolicy& miss, SnpRec& out);
  void ensure_index_warn();
};

} // namespace eqtl
