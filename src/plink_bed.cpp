/* eqtl — PLINK bed: sequential block fread + 2-bit lookup + buffer reuse */
#include "eqtl/plink_bed.hpp"
#include "eqtl/util.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <fstream>

namespace eqtl {

PlinkBed::~PlinkBed() {
  if (bed_fp_) {
    std::fclose(bed_fp_);
    bed_fp_ = nullptr;
  }
}

void PlinkBed::init_lut() {
  pair_lut_[0] = 2;
  pair_lut_[1] = -1;
  pair_lut_[2] = 1;
  pair_lut_[3] = 0;
}

void PlinkBed::read_fam(const std::string& path) {
  std::ifstream in(path);
  if (!in) die("cannot open " + path);
  samples_.clear();
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto t = split_ws(line);
    if (t.size() < 2) die("malformed fam line (need FID IID ...): " + path);
    samples_.push_back(t[1]);
  }
  n_file_ = samples_.size();
  if (n_file_ == 0) die("empty fam: " + path);
}

void PlinkBed::read_bim(const std::string& path) {
  std::ifstream in(path);
  if (!in) die("cannot open " + path);
  sites_.clear();
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto t = split_ws(line);
    if (t.size() < 6) die("malformed bim line: " + path);
    BimSite s;
    s.chrom = t[0];
    s.id = t[1];
    s.cm = std::stod(t[2]);
    s.pos = std::stoll(t[3]);
    s.a1 = t[4];
    s.a2 = t[5];
    sites_.push_back(std::move(s));
  }
  if (sites_.empty()) die("empty bim: " + path);
}

void PlinkBed::build_chrom_ranges() {
  chrom_range_.clear();
  if (sites_.empty()) return;
  size_t i = 0;
  while (i < sites_.size()) {
    const std::string key = chrom_key(sites_[i].chrom);
    size_t j = i + 1;
    while (j < sites_.size() && chrom_equal(sites_[j].chrom, sites_[i].chrom)) ++j;
    chrom_range_[key] = {i, j};
    chrom_range_[sites_[i].chrom] = {i, j};
    i = j;
  }
}

void PlinkBed::open_bed(const std::string& path) {
  if (bed_fp_) {
    std::fclose(bed_fp_);
    bed_fp_ = nullptr;
  }
  bed_fp_ = std::fopen(path.c_str(), "rb");
  if (!bed_fp_) die("cannot open " + path);
  file_buf_.assign(1 << 20, 0);
  std::setvbuf(bed_fp_, file_buf_.data(), _IOFBF, file_buf_.size());

  unsigned char magic[3];
  if (std::fread(magic, 1, 3, bed_fp_) != 3) die("cannot read bed magic: " + path);
  if (magic[0] != 0x6c || magic[1] != 0x1b) die("not a PLINK1 .bed (bad magic): " + path);
  if (magic[2] != 0x01) die("only SNP-major PLINK .bed is supported: " + path);
  bytes_per_snp_ = (n_file_ + 3) / 4;
  block_buf_.assign(kBlockSnps * bytes_per_snp_, 0);

  if (std::fseek(bed_fp_, 0, SEEK_END) != 0) die("bed seek end failed: " + path);
  const long fsz = std::ftell(bed_fp_);
  if (fsz < 0) die("bed ftell failed: " + path);
  const uint64_t expect = 3ull + bytes_per_snp_ * static_cast<uint64_t>(sites_.size());
  if (static_cast<uint64_t>(fsz) < expect) die("bed size too small for bim/fam: " + path);
  if (std::fseek(bed_fp_, 3, SEEK_SET) != 0) die("bed seek data failed: " + path);
}

void PlinkBed::open(const std::string& prefix) {
  prefix_ = prefix;
  init_lut();
  read_fam(prefix + ".fam");
  read_bim(prefix + ".bim");
  build_chrom_ranges();
  open_bed(prefix + ".bed");
  info("bfile: " + std::to_string(n_file_) + " samples, " +
       std::to_string(sites_.size()) + " SNPs [" + prefix + "]");
}

void PlinkBed::set_sample_order(const std::vector<std::string>& sample_ids) {
  const auto m = index_map(samples_);
  sample_col_.clear();
  sample_col_.reserve(sample_ids.size());
  for (const auto& id : sample_ids) {
    auto it = m.find(id);
    if (it == m.end()) die("sample not in fam: " + id);
    sample_col_.push_back(it->second);
  }
  // Pre-size reused dosage buffer
  snp_reuse_.dosage.resize(sample_col_.size());
}

bool PlinkBed::seek_snp(size_t snp_idx) {
  if (!bed_fp_) return false;
  const long off = static_cast<long>(3 + snp_idx * bytes_per_snp_);
  return std::fseek(bed_fp_, off, SEEK_SET) == 0;
}

bool PlinkBed::decode_row(size_t snp_idx, const uint8_t* row, const MissPolicy& miss, double maf_min,
                          SnpRec& out) {
  if (sample_col_.empty() || snp_idx >= sites_.size()) return false;
  const int n_an = static_cast<int>(sample_col_.size());
  if (static_cast<int>(out.dosage.size()) != n_an) out.dosage.resize(static_cast<size_t>(n_an));
  // Mark all missing first via NaN
  std::fill(out.dosage.begin(), out.dosage.end(), std::numeric_limits<double>::quiet_NaN());
  int n_miss = 0;
  double sum = 0.0;
  int n_ok = 0;
  for (int i = 0; i < n_an; ++i) {
    const int col = sample_col_[static_cast<size_t>(i)];
    if (col < 0 || static_cast<size_t>(col) >= n_file_) return false;
    const size_t byte_i = static_cast<size_t>(col) / 4;
    const int j = col % 4;
    const uint8_t b = row[byte_i];
    const int pair = (b >> (2 * j)) & 3;
    const int a1c = pair_lut_[pair];
    if (a1c < 0) {
      ++n_miss;
      continue;
    }
    out.dosage[static_cast<size_t>(i)] = static_cast<double>(a1c);
    sum += a1c;
    ++n_ok;
  }
  if (n_ok == 0) return false;

  const double miss_frac = static_cast<double>(n_miss) / static_cast<double>(n_an);
  if (miss_frac > miss.max_miss + 1e-15) return false;
  if (miss.hand == MissHand::Filter && miss.max_miss <= 0.0 && n_miss > 0) return false;
  if (n_miss > 0) {
    const double mu = sum / n_ok;
    for (int i = 0; i < n_an; ++i)
      if (!std::isfinite(out.dosage[static_cast<size_t>(i)])) out.dosage[static_cast<size_t>(i)] = mu;
  }

  double maf = (sum / static_cast<double>(n_ok)) / 2.0;
  if (maf > 0.5) maf = 1.0 - maf;
  if (maf < 1e-12) return false;
  if (!pass_maf(maf, maf_min)) return false;

  const auto& st = sites_[snp_idx];
  out.chrom = st.chrom;
  out.pos = st.pos;
  out.ref = st.a2;
  out.alt = st.a1;
  out.id = st.id;
  out.maf = maf;
  return true;
}

bool PlinkBed::for_each_range(size_t lo, size_t hi, const MissPolicy& miss, double maf_min,
                              const std::function<bool(const SnpRec&)>& fn) {
  return for_each_range_pos(lo, hi, 0, std::numeric_limits<int64_t>::max(), miss, maf_min, fn);
}

bool PlinkBed::for_each_range_pos(size_t lo, size_t hi, int64_t p0, int64_t p1,
                                  const MissPolicy& miss, double maf_min,
                                  const std::function<bool(const SnpRec&)>& fn) {
  if (lo >= hi || lo >= sites_.size()) return true;
  if (hi > sites_.size()) hi = sites_.size();
  if (!seek_snp(lo)) die("bed seek failed");

  size_t t = lo;
  while (t < hi) {
    const size_t n_take = std::min(kBlockSnps, hi - t);
    const size_t nbytes = n_take * bytes_per_snp_;
    if (std::fread(block_buf_.data(), 1, nbytes, bed_fp_) != nbytes) {
      die("bed fread short at SNP " + std::to_string(t));
    }
    for (size_t k = 0; k < n_take; ++k) {
      const size_t idx = t + k;
      const int64_t pos = sites_[idx].pos;
      if (pos < p0 || pos > p1) continue;
      const uint8_t* row = block_buf_.data() + k * bytes_per_snp_;
      if (!decode_row(idx, row, miss, maf_min, snp_reuse_)) continue;
      if (!fn(snp_reuse_)) return false;
    }
    t += n_take;
  }
  return true;
}

void PlinkBed::for_each_snp(const MissPolicy& miss, double maf_min,
                            const std::function<bool(const SnpRec&)>& fn) {
  for_each_range(0, sites_.size(), miss, maf_min, fn);
}

void PlinkBed::for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                                   const MissPolicy& miss, double maf_min,
                                   const std::function<bool(const SnpRec&)>& fn) {
  if (start < 1) start = 1;
  if (end < start) return;

  size_t clo = 0, chi = 0;
  bool found = false;
  auto it = chrom_range_.find(chrom);
  if (it == chrom_range_.end()) it = chrom_range_.find(chrom_key(chrom));
  if (it != chrom_range_.end()) {
    clo = it->second.first;
    chi = it->second.second;
    found = true;
  } else {
    for (const auto& kv : chrom_range_) {
      if (kv.second.first < sites_.size() && chrom_equal(sites_[kv.second.first].chrom, chrom)) {
        clo = kv.second.first;
        chi = kv.second.second;
        found = true;
        break;
      }
    }
  }
  if (!found) return;

  auto pos_at = [&](size_t i) { return sites_[i].pos; };
  size_t lo = clo, hi = chi;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (pos_at(mid) < start) lo = mid + 1;
    else hi = mid;
  }
  size_t left = lo;
  lo = clo;
  hi = chi;
  while (lo < hi) {
    size_t mid = lo + (hi - lo) / 2;
    if (pos_at(mid) <= end) lo = mid + 1;
    else hi = mid;
  }
  size_t right = lo;
  if (left >= right) return;
  for_each_range_pos(left, right, start, end, miss, maf_min, fn);
}

std::vector<SnpRec> PlinkBed::load_all(const MissPolicy& miss, double maf_min, int max_snps) {
  std::vector<SnpRec> all;
  for_each_snp(miss, maf_min, [&](const SnpRec& s) {
    all.push_back(s);
    return !(max_snps > 0 && static_cast<int>(all.size()) >= max_snps);
  });
  return all;
}

std::vector<SnpRec> PlinkBed::load_region(const std::string& chrom, int64_t start, int64_t end,
                                          const MissPolicy& miss, double maf_min) {
  std::vector<SnpRec> out;
  out.reserve(4096);
  for_each_snp_region(chrom, start, end, miss, maf_min, [&](const SnpRec& s) {
    out.push_back(s);
    return true;
  });
  return out;
}

} // namespace eqtl
