#define _FILE_OFFSET_BITS 64
/* eqtl — PLINK bed: sequential block fread + 2-bit lookup + buffer reuse */
#include "eqtl/plink_bed.hpp"
#include <unordered_set>
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
  if (samples_.empty()) die("empty fam: " + path);
  assert_unique_ids(samples_, "plink fam IID");
  n_file_ = samples_.size();
}

void PlinkBed::read_bim(const std::string& path) {
  FILE* fp = std::fopen(path.c_str(), "rb");
  if (!fp) die("cannot open " + path);

  sites_.clear();
  sites_.reserve(1 << 20);

  // 4MB stdio buffer for fast streaming I/O
  std::vector<char> file_buf(4 << 20);
  std::setvbuf(fp, file_buf.data(), _IOFBF, file_buf.size());

  // Line buffer for streaming
  std::vector<char> line_buf(1024);
  while (true) {
    // Read one line
    if (std::fgets(line_buf.data(), static_cast<int>(line_buf.size()), fp) == nullptr) break;
    // Grow line buffer if line was truncated (no newline at end)
    while (line_buf.back() != '\0' && line_buf.back() != '\n' && !std::feof(fp)) {
      size_t old_sz = line_buf.size();
      line_buf.resize(old_sz * 2);
      if (std::fgets(line_buf.data() + old_sz, static_cast<int>(old_sz + 1), fp) == nullptr) break;
    }

    const char* p = line_buf.data();
    // Skip leading whitespace
    while (*p == ' ' || *p == '\t') ++p;
    // Skip empty lines
    if (*p == '\n' || *p == '\r' || *p == '\0') continue;

    // Parse 6 fields: chrom id cm pos a1 a2
    auto read_field = [&]() -> std::string {
      const char* start = p;
      while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') ++p;
      return std::string(start, p);
    };
    auto skip_ws = [&]() { while (*p == ' ' || *p == '\t') ++p; };
    auto parse_int64_checked = [&](const char* field_name) -> int64_t {
      const char* start = p;
      int64_t v = 0;
      bool neg = false;
      if (*p == '-') { neg = true; ++p; }
      if (*p < '0' || *p > '9') die("malformed bim " + std::string(field_name) + " in: " + path);
      while (*p >= '0' && *p <= '9') { v = v * 10 + (*p - '0'); ++p; }
      if (p == start || (p == start + 1 && neg)) die("malformed bim " + std::string(field_name) + " in: " + path);
      return neg ? -v : v;
    };
    auto parse_double_checked = [&](const char* field_name) -> double {
      const char* start = p;
      double v = 0;
      bool neg = false;
      if (*p == '-') { neg = true; ++p; }
      if (*p < '0' || *p > '9') die("malformed bim " + std::string(field_name) + " in: " + path);
      while (*p >= '0' && *p <= '9') { v = v * 10 + (*p - '0'); ++p; }
      if (*p == '.') {
        ++p;
        double frac = 1;
        while (*p >= '0' && *p <= '9') { frac *= 0.1; v += (*p - '0') * frac; ++p; }
      }
      if (p == start || (p == start + 1 && neg)) die("malformed bim " + std::string(field_name) + " in: " + path);
      return neg ? -v : v;
    };

    BimSite s;
    s.chrom = read_field();
    skip_ws();
    s.id = read_field();
    skip_ws();
    s.cm = parse_double_checked("cm");
    skip_ws();
    s.pos = parse_int64_checked("pos");
    skip_ws();
    s.a1 = read_field();
    skip_ws();
    s.a2 = read_field();

    if (s.chrom.empty() || s.id.empty() || s.a1.empty() || s.a2.empty())
      die("malformed bim line in: " + path);
    sites_.push_back(std::move(s));
  }
  std::fclose(fp);

  if (sites_.empty()) die("empty bim: " + path);
  // ponytail: string key once at load; die on identical chrom:pos:A1:A2
  {
    std::unordered_set<std::string> seen;
    seen.reserve(sites_.size() * 2);
    for (const auto& s : sites_) {
      std::string k;
      k.reserve(s.chrom.size() + s.a1.size() + s.a2.size() + 24);
      k.append(s.chrom).push_back('\t');
      k.append(std::to_string(s.pos)).push_back('\t');
      k.append(s.a1).push_back('\t');
      k.append(s.a2);
      if (!seen.insert(std::move(k)).second)
        die("plink bim site: duplicate id: " + s.chrom + ":" + std::to_string(s.pos) + ":" + s.a1 +
            ":" + s.a2);
    }
  }
  // for_each_snp_region binary-searches positions inside one contiguous block
  // per chromosome; an interleaved or unsorted bim silently drops SNPs from
  // cis windows (verified: shuffled bim loses ~37% of cis pairs, no warning).
  {
    std::unordered_set<std::string> seen_chroms;
    for (size_t i = 0; i < sites_.size(); ++i) {
      const bool first_of_run =
          i == 0 || !chrom_equal(sites_[i - 1].chrom, sites_[i].chrom);
      if (first_of_run) {
        if (!seen_chroms.insert(chrom_key(sites_[i].chrom)).second)
          die("plink bim: chromosome " + sites_[i].chrom + " appears in more than one "
              "block (line " + std::to_string(i + 1) + " of " + path +
              "); region queries would miss SNPs — re-sort: plink --bfile <prefix> --make-bed");
      } else if (sites_[i - 1].pos > sites_[i].pos) {
        die("plink bim: positions not sorted within chromosome " + sites_[i].chrom +
            " (line " + std::to_string(i + 1) + " of " + path +
            "); region queries would miss SNPs — re-sort: plink --bfile <prefix> --make-bed");
      }
    }
  }
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

std::vector<std::string> PlinkBed::chromosomes() const {
  std::vector<std::string> out;
  if (sites_.empty()) return out;
  out.push_back(sites_[0].chrom);
  for (size_t i = 1; i < sites_.size(); ++i) {
    if (!chrom_equal(sites_[i].chrom, out.back()))
      out.push_back(sites_[i].chrom);
  }
  return out;
}

void PlinkBed::open_bed(const std::string& path) {
  if (bed_fp_) {
    std::fclose(bed_fp_);
    bed_fp_ = nullptr;
  }
  bed_fp_ = std::fopen(path.c_str(), "rb");
  if (!bed_fp_) die("cannot open " + path);
  file_buf_.assign(4 << 20, 0); // 4MB stdio buffer
  std::setvbuf(bed_fp_, file_buf_.data(), _IOFBF, file_buf_.size());

  unsigned char magic[3];
  if (std::fread(magic, 1, 3, bed_fp_) != 3) die("cannot read bed magic: " + path);
  if (magic[0] != 0x6c || magic[1] != 0x1b) die("not a PLINK1 .bed (bad magic): " + path);
  if (magic[2] != 0x01) die("only SNP-major PLINK .bed is supported: " + path);
  bytes_per_snp_ = (n_file_ + 3) / 4;
  // ponytail: cap block by ~4MB; floor 4096 SNPs
  block_snps_ = std::max(kMinBlockSnps, kTargetBlockBytes / std::max<size_t>(bytes_per_snp_, 1));
  block_buf_.assign(block_snps_ * bytes_per_snp_, 0);

  if (fseeko(bed_fp_, 0, SEEK_END) != 0) die("bed seek end failed: " + path);
  const off_t fsz = ftello(bed_fp_);
  if (fsz < 0) die("bed ftell failed: " + path);
  const uint64_t expect = 3ull + bytes_per_snp_ * static_cast<uint64_t>(sites_.size());
  if (static_cast<uint64_t>(fsz) < expect) die("bed size too small for bim/fam: " + path);
  if (fseeko(bed_fp_, 3, SEEK_SET) != 0) die("bed seek data failed: " + path);
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
  snp_reuse_.dosage.resize(sample_col_.size());
}

bool PlinkBed::seek_snp(size_t snp_idx) {
  if (!bed_fp_) return false;
  const off_t off = static_cast<off_t>(3 + snp_idx * bytes_per_snp_);
  return fseeko(bed_fp_, off, SEEK_SET) == 0;
}

bool PlinkBed::decode_row(size_t snp_idx, const uint8_t* row, const MissPolicy& miss, double maf_min,
                          SnpRec& out) {
  if (sample_col_.empty() || snp_idx >= sites_.size()) return false;
  const int n_an = static_cast<int>(sample_col_.size());
  if (static_cast<int>(out.dosage.size()) != n_an) out.dosage.resize(static_cast<size_t>(n_an));
  int n_miss = 0;
  double sum = 0.0;
  int n_ok = 0;
  // ponytail: write dosage directly; miss positions set after if needed (skip full NaN memset)
  for (int i = 0; i < n_an; ++i) {
    const int col = sample_col_[static_cast<size_t>(i)];
    if (col < 0 || static_cast<size_t>(col) >= n_file_) return false;
    const size_t byte_i = static_cast<size_t>(col) / 4;
    const int j = col % 4;
    const uint8_t b = row[byte_i];
    const int pair = (b >> (2 * j)) & 3;
    const int a1c = pair_lut_[pair];
    if (a1c < 0) {
      out.dosage[static_cast<size_t>(i)] = std::numeric_limits<double>::quiet_NaN();
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
  // filter: any remaining missing → drop SNP (complete-case; no silent NaN into LM)
  if (miss.hand == MissHand::Filter && n_miss > 0) return false;
  if (n_miss > 0 && miss.hand == MissHand::Impute) {
    const double mu = sum / n_ok;
    for (int i = 0; i < n_an; ++i)
      if (!std::isfinite(out.dosage[static_cast<size_t>(i)])) out.dosage[static_cast<size_t>(i)] = mu;
  }

  // effect-allele AF = mean(A1 count)/2 (diploid bed). Not folded to minor.
  double af = (sum / static_cast<double>(n_ok)) / 2.0;
  if (af < 0) af = 0;
  if (af > 1) af = 1;
  if (af < 1e-12 || af > 1.0 - 1e-12) return false;
  if (!pass_maf(af, maf_min)) return false; // pass_maf uses min(af,1-af) internally

  const auto& st = sites_[snp_idx];
  out.chrom = st.chrom;
  out.pos = st.pos;
  // effect allele = A1 → alt; other allele A2 → ref (PLINK convention, not genome REF/ALT)
  out.ref = st.a2;
  out.alt = st.a1;
  out.id = st.id;
  out.maf = af;
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
    const size_t n_take = std::min(block_snps_, hi - t);
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
