/* eqtl — PLINK bed/bim/fam reader (encoding matches GEMMA ReadFile_bed) */
#include "eqtl/plink_bed.hpp"
#include "eqtl/util.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

namespace eqtl {

PlinkBed::~PlinkBed() {
  if (bed_.is_open()) bed_.close();
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
    samples_.push_back(t[1]); // IID
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

void PlinkBed::open_bed(const std::string& path) {
  bed_.open(path, std::ios::binary);
  if (!bed_) die("cannot open " + path);
  unsigned char magic[3];
  bed_.read(reinterpret_cast<char*>(magic), 3);
  if (!bed_ || magic[0] != 0x6c || magic[1] != 0x1b) {
    die("not a PLINK1 .bed (bad magic): " + path);
  }
  // 0x01 = SNP-major (required); 0x00 = individual-major
  snp_major_ = (magic[2] == 0x01);
  if (!snp_major_) die("only SNP-major PLINK .bed is supported: " + path);
  bytes_per_snp_ = (n_file_ + 3) / 4;
  rowbuf_.assign(bytes_per_snp_, 0);
  // verify size roughly
  bed_.seekg(0, std::ios::end);
  const auto fsz = static_cast<uint64_t>(bed_.tellg());
  const uint64_t expect = 3ull + bytes_per_snp_ * static_cast<uint64_t>(sites_.size());
  if (fsz < expect) {
    die("bed size too small for bim/fam: " + path);
  }
  bed_.clear();
  bed_.seekg(3, std::ios::beg);
}

void PlinkBed::open(const std::string& prefix) {
  prefix_ = prefix;
  read_fam(prefix + ".fam");
  read_bim(prefix + ".bim");
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
}

// PLINK SNP-major byte: 4 genotypes, low pair first.
// GEMMA: 00→2 (A1/A1), 10→1 (het), 11→0 (A2/A2), 01→miss  [bitset bit order]
// Equivalent: bits (b0,b1) as integer b0+2*b1: 0→2, 2→1, 3→0, 1→miss
static inline int plink_pair_to_a1_count(int pair) {
  // pair = bit0 | (bit1<<1)
  switch (pair) {
    case 0: return 2; // 00 A1/A1
    case 2: return 1; // 10 het
    case 3: return 0; // 11 A2/A2
    case 1: return -1; // 01 missing
    default: return -1;
  }
}

bool PlinkBed::decode_snp(size_t snp_idx, const MissPolicy& miss, double maf_min, SnpRec& out) {
  if (sample_col_.empty()) return false;
  bed_.seekg(static_cast<std::streamoff>(3 + snp_idx * bytes_per_snp_), std::ios::beg);
  bed_.read(reinterpret_cast<char*>(rowbuf_.data()),
            static_cast<std::streamsize>(bytes_per_snp_));
  if (!bed_ || static_cast<size_t>(bed_.gcount()) != bytes_per_snp_) return false;

  const int n_an = static_cast<int>(sample_col_.size());
  out.dosage.assign(n_an, std::numeric_limits<double>::quiet_NaN());
  int n_miss = 0;
  double sum = 0.0;
  int n_ok = 0;
  for (int i = 0; i < n_an; ++i) {
    const int col = sample_col_[i];
    if (col < 0 || static_cast<size_t>(col) >= n_file_) return false;
    const size_t byte_i = static_cast<size_t>(col) / 4;
    const int j = col % 4; // genotype slot in byte
    const uint8_t b = rowbuf_[byte_i];
    const int pair = (b >> (2 * j)) & 3;
    const int a1c = plink_pair_to_a1_count(pair);
    if (a1c < 0) {
      ++n_miss;
      continue;
    }
    out.dosage[i] = static_cast<double>(a1c);
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
      if (!std::isfinite(out.dosage[i])) out.dosage[i] = mu;
  }

  // MAF on non-missing before impute would use sum/n_ok; after impute uses full N.
  // GEMMA uses non-missing: maf = sum/(2*(ni_test-n_miss)). Match that using sum/n_ok.
  double maf = (sum / static_cast<double>(n_ok)) / 2.0;
  if (maf > 0.5) maf = 1.0 - maf;
  if (maf < 1e-12) return false;
  if (!pass_maf(maf, maf_min)) return false;

  const auto& st = sites_[snp_idx];
  out.chrom = st.chrom;
  out.pos = st.pos;
  out.ref = st.a2; // A2 as "other"; dosage is A1 count (effect allele A1)
  out.alt = st.a1;
  out.id = st.id;
  out.maf = maf;
  return true;
}

void PlinkBed::for_each_snp(const MissPolicy& miss, double maf_min,
                            const std::function<bool(const SnpRec&)>& fn) {
  SnpRec snp;
  for (size_t t = 0; t < sites_.size(); ++t) {
    if (!decode_snp(t, miss, maf_min, snp)) continue;
    if (!fn(snp)) break;
  }
}

void PlinkBed::for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                                   const MissPolicy& miss, double maf_min,
                                   const std::function<bool(const SnpRec&)>& fn) {
  if (start < 1) start = 1;
  if (end < start) return;
  SnpRec snp;
  for (size_t t = 0; t < sites_.size(); ++t) {
    const auto& st = sites_[t];
    if (!chrom_equal(st.chrom, chrom) || st.pos < start || st.pos > end) continue;
    if (!decode_snp(t, miss, maf_min, snp)) continue;
    if (!fn(snp)) break;
  }
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
