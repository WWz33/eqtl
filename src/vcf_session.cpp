/* eqtl — VCF/BCF session, GT dosage */
#include "eqtl/vcf_session.hpp"
#include "eqtl/util.hpp"
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <cmath>
#include <limits>
#include <cstdlib>

namespace eqtl {

VcfSession::~VcfSession() {
  if (tbx_) tbx_destroy(static_cast<tbx_t*>(tbx_));
  if (idx_) hts_idx_destroy(static_cast<hts_idx_t*>(idx_));
  if (hdr_) bcf_hdr_destroy(static_cast<bcf_hdr_t*>(hdr_));
  if (fp_) hts_close(static_cast<htsFile*>(fp_));
  tbx_ = idx_ = hdr_ = fp_ = nullptr;
}

void VcfSession::open(const std::string& path) {
  path_ = path;
  fp_ = hts_open(path.c_str(), "r");
  if (!fp_) die("cannot open VCF/BCF: " + path);
  hdr_ = bcf_hdr_read(static_cast<htsFile*>(fp_));
  if (!hdr_) die("cannot read VCF header: " + path);
  auto* hdr = static_cast<bcf_hdr_t*>(hdr_);
  const int ns = bcf_hdr_nsamples(hdr);
  samples_.clear();
  samples_.reserve(ns);
  for (int i = 0; i < ns; ++i) samples_.emplace_back(hdr->samples[i]);
  contigs_.clear();
  for (int i = 0; i < hdr->n[BCF_DT_CTG]; ++i)
    contigs_.emplace_back(hdr->id[BCF_DT_CTG][i].key);

  idx_ = bcf_index_load(path.c_str());
  if (!idx_) tbx_ = tbx_index_load(path.c_str());
  indexed_ = (idx_ != nullptr) || (tbx_ != nullptr);
  if (!indexed_) ensure_index_warn();
  info("vcf: " + std::to_string(samples_.size()) + " samples, " +
       std::to_string(contigs_.size()) + " contigs" +
       (indexed_ ? " (indexed)" : " (no index)"));
}

void VcfSession::ensure_index_warn() {
  if (warned_no_index_) return;
  warned_no_index_ = true;
  warn("VCF has no CSI/TBI; sequential scan is slower. Run: bcftools index -t <vcf>");
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
}

bool VcfSession::parse_record(void* rec_v, MissHand miss, SnpRec& out) {
  auto* rec = static_cast<bcf1_t*>(rec_v);
  // Must use the header that read this record (active_hdr_ or hdr_).
  auto* hdr = static_cast<bcf_hdr_t*>(active_hdr_ ? active_hdr_ : hdr_);
  if (bcf_unpack(rec, BCF_UN_STR | BCF_UN_FMT) < 0) return false;
  if (rec->n_allele != 2) return false;

  int ngt = 0;
  int32_t* gt = nullptr;
  if (bcf_get_genotypes(hdr, rec, &gt, &ngt) <= 0 || !gt) {
    if (gt) free(gt);
    return false;
  }
  const int n_an = static_cast<int>(sample_col_.size());
  if (n_an == 0) {
    free(gt);
    return false;
  }
  const int ns_all = bcf_hdr_nsamples(hdr);
  if (ns_all <= 0 || ngt % ns_all != 0) {
    free(gt);
    return false;
  }
  const int max_pl = ngt / ns_all;
  out.dosage.assign(n_an, std::numeric_limits<double>::quiet_NaN());
  int n_miss = 0;
  double sum = 0.0;
  int n_ok = 0;
  for (int i = 0; i < n_an; ++i) {
    const int col = sample_col_[i];
    if (col < 0 || col >= ns_all) {
      free(gt);
      return false;
    }
    const int32_t a0 = gt[col * max_pl];
    const int32_t a1 = (max_pl > 1) ? gt[col * max_pl + 1] : bcf_int32_vector_end;
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
    out.dosage[i] = static_cast<double>(d);
    sum += d;
    ++n_ok;
  }
  free(gt);
  if (n_ok == 0) return false;
  if (miss == MissHand::Filter && n_miss > 0) return false;
  if (miss == MissHand::Impute && n_miss > 0) {
    const double mu = sum / n_ok;
    for (int i = 0; i < n_an; ++i)
      if (!std::isfinite(out.dosage[i])) out.dosage[i] = mu;
  }
  double mean = 0.0;
  for (double d : out.dosage) mean += d;
  mean /= n_an;
  out.maf = mean / 2.0;
  if (out.maf > 0.5) out.maf = 1.0 - out.maf;
  if (out.maf < 1e-12) return false;
  out.chrom = bcf_hdr_id2name(hdr, rec->rid);
  out.pos = rec->pos + 1;
  out.ref = rec->d.allele[0] ? rec->d.allele[0] : ".";
  out.alt = rec->d.allele[1] ? rec->d.allele[1] : ".";
  if (rec->d.id && rec->d.id[0] && std::string(rec->d.id) != ".")
    out.id = rec->d.id;
  else
    out.id = out.chrom + ":" + std::to_string(out.pos) + ":" + out.ref + ":" + out.alt;
  return true;
}

void VcfSession::for_each_snp(MissHand miss, const std::function<bool(const SnpRec&)>& fn) {
  // Use original fp/hdr: seek to first record after header via reopen is safer;
  // set active_hdr_ to the header used for bcf_read.
  htsFile* rfp = hts_open(path_.c_str(), "r");
  if (!rfp) die("reopen VCF failed: " + path_);
  bcf_hdr_t* rh = bcf_hdr_read(rfp);
  active_hdr_ = rh;
  bcf1_t* rec = bcf_init();
  SnpRec snp;
  while (bcf_read(rfp, rh, rec) == 0) {
    if (!parse_record(rec, miss, snp)) continue;
    if (!fn(snp)) break;
  }
  bcf_destroy(rec);
  active_hdr_ = nullptr;
  bcf_hdr_destroy(rh);
  hts_close(rfp);
}

void VcfSession::for_each_snp_region(const std::string& chrom, int64_t start, int64_t end,
                                     MissHand miss, const std::function<bool(const SnpRec&)>& fn) {
  for_each_snp(miss, [&](const SnpRec& s) {
    if (s.chrom != chrom || s.pos < start || s.pos > end) return true;
    return fn(s);
  });
}

std::vector<SnpRec> VcfSession::load_all(MissHand miss, int max_snps) {
  std::vector<SnpRec> all;
  for_each_snp(miss, [&](const SnpRec& s) {
    all.push_back(s);
    return !(max_snps > 0 && static_cast<int>(all.size()) >= max_snps);
  });
  return all;
}

}  // namespace eqtl
