/* eqtl — VCF/BCF session, GT dosage (shared htsFile; no per-query reopen) */
#include "eqtl/vcf_session.hpp"
#include "eqtl/util.hpp"
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/bgzf.h>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace eqtl {

VcfSession::~VcfSession() {
  if (gt_) free(gt_);
  if (rec_) bcf_destroy(static_cast<bcf1_t*>(rec_));
  if (tbx_) tbx_destroy(static_cast<tbx_t*>(tbx_));
  if (idx_) hts_idx_destroy(static_cast<hts_idx_t*>(idx_));
  if (hdr_) bcf_hdr_destroy(static_cast<bcf_hdr_t*>(hdr_));
  if (fp_) hts_close(static_cast<htsFile*>(fp_));
  gt_ = nullptr;
  ngt_ = 0;
  rec_ = tbx_ = idx_ = hdr_ = fp_ = nullptr;
}

void VcfSession::open(const std::string& path) {
  path_ = path;
  fp_ = hts_open(path.c_str(), "r");
  if (!fp_) die("cannot open VCF/BCF: " + path);
  hdr_ = bcf_hdr_read(static_cast<htsFile*>(fp_));
  if (!hdr_) die("cannot read VCF header: " + path);

  // Position after header for full-scan rewind (bgzf / plain).
  data_off_ = 0;
  BGZF* bgz = hts_get_bgzfp(static_cast<htsFile*>(fp_));
  if (bgz) {
    data_off_ = bgzf_tell(bgz);
  }

  auto* hdr = static_cast<bcf_hdr_t*>(hdr_);
  const int ns = bcf_hdr_nsamples(hdr);
  samples_.clear();
  samples_.reserve(ns);
  for (int i = 0; i < ns; ++i) samples_.emplace_back(hdr->samples[i]);
  contigs_.clear();
  for (int i = 0; i < hdr->n[BCF_DT_CTG]; ++i)
    contigs_.emplace_back(hdr->id[BCF_DT_CTG][i].key);

  // Index: TBI for VCF.gz, CSI for BCF (bcftools index). Order reduces htslib noise.
  const bool likely_vcf = path.size() >= 4 &&
      (path.rfind(".vcf") != std::string::npos || path.rfind(".VCF") != std::string::npos);
  if (likely_vcf) {
    tbx_ = tbx_index_load(path.c_str());
    if (!tbx_) {
      hts_idx_t* csi = bcf_index_load(path.c_str());
      if (csi) idx_ = csi;
    }
  } else {
    hts_idx_t* csi = bcf_index_load(path.c_str());
    if (csi) {
      idx_ = csi;
    } else {
      tbx_ = tbx_index_load(path.c_str());
    }
  }
  indexed_ = (idx_ != nullptr) || (tbx_ != nullptr);
  if (!indexed_) ensure_index_warn();

  rec_ = bcf_init();
  if (!rec_) die("bcf_init failed");

  info("vcf: " + std::to_string(samples_.size()) + " samples, " +
       std::to_string(contigs_.size()) + " contigs" +
       (indexed_ ? " (indexed)" : " (no index)"));
}

void VcfSession::ensure_index_warn() {
  if (warned_no_index_) return;
  warned_no_index_ = true;
  warn("VCF/BCF has no CSI/TBI; region queries scan the whole file. "
       "Index with: bcftools index -t <file.vcf.gz>  or  bcftools index <file.bcf>");
}

bool VcfSession::rewind_to_data() {
  auto* fp = static_cast<htsFile*>(fp_);
  if (!fp) return false;
  BGZF* bgz = hts_get_bgzfp(fp);
  if (bgz) {
    return bgzf_seek(bgz, data_off_, SEEK_SET) == 0;
  }
  // Uncompressed / non-bgzf: re-open and discard a fresh header read (keep hdr_).
  hts_close(fp);
  fp_ = hts_open(path_.c_str(), "r");
  if (!fp_) return false;
  bcf_hdr_t* tmp = bcf_hdr_read(static_cast<htsFile*>(fp_));
  if (!tmp) return false;
  bcf_hdr_destroy(tmp);
  return true;
}

void VcfSession::set_sample_order(const std::vector<std::string>& sample_ids) {
  const auto m = index_map(samples_);
  sample_col_.clear();
  sample_col_.reserve(sample_ids.size());
  for (const auto& id : sample_ids) {
    const auto it = m.find(id);
    if (it == m.end()) die("sample not in VCF: " + id);
    sample_col_.push_back(it->second);
  }
  snp_reuse_.dosage.resize(sample_col_.size());
}

std::string VcfSession::resolve_contig(const std::string& chrom) const {
  for (const auto& c : contigs_) {
    if (chrom_equal(c, chrom)) return c;
  }
  return {};
}

bool VcfSession::parse_record(void* rec_v, const MissPolicy& miss, SnpRec& out) {
  auto* rec = static_cast<bcf1_t*>(rec_v);
  auto* hdr = static_cast<bcf_hdr_t*>(hdr_);
  if (bcf_unpack(rec, BCF_UN_STR | BCF_UN_FMT) < 0) return false;
  if (rec->n_allele != 2) return false;

  // Reuse gt_ buffer (bcf_get_genotypes reallocs as needed)
  if (bcf_get_genotypes(hdr, rec, &gt_, &ngt_) <= 0 || !gt_) {
    return false;
  }
  const int n_an = static_cast<int>(sample_col_.size());
  if (n_an == 0) return false;
  const int ns_all = bcf_hdr_nsamples(hdr);
  if (ns_all <= 0 || ngt_ % ns_all != 0) return false;
  const int max_pl = ngt_ / ns_all;
  if (static_cast<int>(out.dosage.size()) != n_an) out.dosage.resize(static_cast<size_t>(n_an));
  std::fill(out.dosage.begin(), out.dosage.end(), std::numeric_limits<double>::quiet_NaN());
  int n_miss = 0;
  double sum = 0.0;
  int n_ok = 0;
  for (int i = 0; i < n_an; ++i) {
    const int col = sample_col_[i];
    if (col < 0 || col >= ns_all) return false;
    const int32_t a0 = gt_[col * max_pl];
    const int32_t a1 = (max_pl > 1) ? gt_[col * max_pl + 1] : bcf_int32_vector_end;
    bool miss_gt = bcf_gt_is_missing(a0);
    if (!miss_gt && max_pl > 1 && a1 != bcf_int32_vector_end && a1 != bcf_int32_missing)
      miss_gt = bcf_gt_is_missing(a1);
    if (miss_gt) {
      ++n_miss;
      continue;
    }
    int d = 0;
    if (!bcf_gt_is_missing(a0) && a0 != bcf_int32_vector_end) d += bcf_gt_allele(a0);
    if (max_pl > 1 && a1 != bcf_int32_vector_end && a1 != bcf_int32_missing && !bcf_gt_is_missing(a1))
      d += bcf_gt_allele(a1);
    out.dosage[static_cast<size_t>(i)] = static_cast<double>(d);
    sum += d;
    ++n_ok;
  }
  if (n_ok == 0) return false;

  const double miss_frac = static_cast<double>(n_miss) / static_cast<double>(n_an);
  if (miss_frac > miss.max_miss + 1e-15) return false;
  if (miss.hand == MissHand::Filter && miss.max_miss <= 0.0 && n_miss > 0) return false;

  double maf = (sum / static_cast<double>(n_ok)) / 2.0;
  if (maf > 0.5) maf = 1.0 - maf;
  if (maf < 1e-12) return false;

  if (n_miss > 0) {
    const double mu = sum / n_ok;
    for (int i = 0; i < n_an; ++i)
      if (!std::isfinite(out.dosage[static_cast<size_t>(i)]))
        out.dosage[static_cast<size_t>(i)] = mu;
  }

  out.maf = maf;
  out.chrom = bcf_hdr_id2name(hdr, rec->rid);
  out.pos = rec->pos + 1;
  out.ref = rec->d.allele[0] ? rec->d.allele[0] : ".";
  out.alt = rec->d.allele[1] ? rec->d.allele[1] : ".";
  if (rec->d.id && rec->d.id[0] && std::strcmp(rec->d.id, ".") != 0)
    out.id = rec->d.id;
  else
    out.id = out.chrom + ":" + std::to_string(out.pos) + ":" + out.ref + ":" + out.alt;
  return true;
}

void VcfSession::for_each_snp(const MissPolicy& miss, double maf_min,
                              const std::function<bool(const SnpRec&)>& fn) {
  if (!fp_ || !hdr_ || !rec_) die("VCF not open");
  if (!rewind_to_data()) die("cannot rewind VCF for full scan: " + path_);
  auto* rfp = static_cast<htsFile*>(fp_);
  auto* rh = static_cast<bcf_hdr_t*>(hdr_);
  auto* rec = static_cast<bcf1_t*>(rec_);
  int rret = 0;
  while ((rret = bcf_read(rfp, rh, rec)) == 0) {
    if (!parse_record(rec, miss, snp_reuse_)) continue;
    if (!pass_maf(snp_reuse_.maf, maf_min)) continue;
    if (!fn(snp_reuse_)) break;
  }
  if (rret < -1) die("VCF read error (bcf_read=" + std::to_string(rret) + "): " + path_);
}

void VcfSession::for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                                     const MissPolicy& miss, double maf_min,
                                     const std::function<bool(const SnpRec&)>& fn) {
  if (start < 1) start = 1;
  if (end < start) return;
  if (!fp_ || !hdr_ || !rec_) die("VCF not open");

  if (!indexed_) {
    ensure_index_warn();
    for_each_snp(miss, maf_min, [&](const SnpRec& s) {
      if (!chrom_equal(s.chrom, chrom) || s.pos < start || s.pos > end) return true;
      return fn(s);
    });
    return;
  }

  const std::string contig = resolve_contig(chrom);
  if (contig.empty()) return;

  const std::string reg = contig + ":" + std::to_string(start) + "-" + std::to_string(end);
  auto* rfp = static_cast<htsFile*>(fp_);
  auto* rh = static_cast<bcf_hdr_t*>(hdr_);
  auto* rec = static_cast<bcf1_t*>(rec_);

  hts_itr_t* itr = nullptr;
  bool use_bcf_itr = false;
  if (idx_) {
    itr = bcf_itr_querys(static_cast<hts_idx_t*>(idx_), rh, reg.c_str());
    use_bcf_itr = (itr != nullptr);
  }
  if (!itr && tbx_) {
    itr = tbx_itr_querys(static_cast<tbx_t*>(tbx_), reg.c_str());
    use_bcf_itr = false;
  }
  if (!itr) {
    for_each_snp(miss, maf_min, [&](const SnpRec& s) {
      if (!chrom_equal(s.chrom, chrom) || s.pos < start || s.pos > end) return true;
      return fn(s);
    });
    return;
  }

  int ret = 0;
  if (idx_ && use_bcf_itr) {
    while ((ret = bcf_itr_next(rfp, itr, rec)) >= 0) {
      if (!parse_record(rec, miss, snp_reuse_)) continue;
      if (!pass_maf(snp_reuse_.maf, maf_min)) continue;
      if (!fn(snp_reuse_)) break;
    }
    if (ret < -1) {
      hts_itr_destroy(itr);
      die("VCF region read error (bcf_itr_next=" + std::to_string(ret) + "): " + path_);
    }
  } else {
    kstring_t sstr = {0, 0, nullptr};
    while ((ret = tbx_itr_next(rfp, static_cast<tbx_t*>(tbx_), itr, &sstr)) >= 0) {
      if (vcf_parse1(&sstr, rh, rec) < 0) continue;
      if (!parse_record(rec, miss, snp_reuse_)) continue;
      if (!pass_maf(snp_reuse_.maf, maf_min)) continue;
      if (!fn(snp_reuse_)) break;
    }
    free(sstr.s);
    if (ret < -1) {
      hts_itr_destroy(itr);
      die("VCF region read error (tbx_itr_next=" + std::to_string(ret) + "): " + path_);
    }
  }
  hts_itr_destroy(itr);
}

std::vector<SnpRec> VcfSession::load_all(const MissPolicy& miss, double maf_min, int max_snps) {
  std::vector<SnpRec> all;
  if (max_snps > 0) all.reserve(static_cast<size_t>(max_snps));
  for_each_snp(miss, maf_min, [&](const SnpRec& s) {
    all.push_back(s);
    return !(max_snps > 0 && static_cast<int>(all.size()) >= max_snps);
  });
  return all;
}

std::vector<SnpRec> VcfSession::load_region(const std::string& chrom, int64_t start, int64_t end,
                                            const MissPolicy& miss, double maf_min) {
  std::vector<SnpRec> out;
  out.reserve(4096);
  for_each_snp_region(chrom, start, end, miss, maf_min, [&](const SnpRec& s) {
    out.push_back(s);
    return true;
  });
  return out;
}

}  // namespace eqtl
