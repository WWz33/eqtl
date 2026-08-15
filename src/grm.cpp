#include "eqtl/grm.hpp"
#include "eqtl/plink_bed.hpp"
#include "eqtl/util.hpp"
#include <fstream>
#include <cstdio>
#include <cmath>

namespace eqtl {

Grm load_grm_gcta(const std::string& prefix) {
  Grm g;
  std::ifstream idf(prefix + ".grm.id");
  if (!idf) die("cannot open " + prefix + ".grm.id");
  std::string a, b;
  // .grm.id: one sample per line; FID IID → use IID, else single token
  while (idf >> a) {
    std::string rest;
    std::getline(idf, rest);
    auto t = split_ws(a + rest);
    if (t.size() >= 2) g.ids.push_back(t[1]); // IID
    else if (t.size() == 1) g.ids.push_back(t[0]);
  }
  const int n = static_cast<int>(g.ids.size());
  if (n == 0) die("empty GRM id file");
  const uint64_t ntri = static_cast<uint64_t>(n) * (n + 1) / 2;
  FILE* fp = std::fopen((prefix + ".grm.bin").c_str(), "rb");
  if (!fp) die("cannot open " + prefix + ".grm.bin");
  std::vector<char> iobuf(1 << 20);
  std::setvbuf(fp, iobuf.data(), _IOFBF, iobuf.size());
  std::vector<float> buf(ntri);
  if (std::fread(buf.data(), sizeof(float), ntri, fp) != ntri)
    die("GRM bin size mismatch");
  std::fclose(fp);
  g.K = Eigen::MatrixXd::Zero(n, n);
  size_t k = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double v = buf[k++];
      g.K(i, j) = v;
      g.K(j, i) = v;
    }
  }
  info("grm: " + std::to_string(n) + " samples from " + prefix);
  return g;
}

void write_grm_gcta(const std::string& prefix, const Grm& g) {
  const int n = static_cast<int>(g.ids.size());
  std::ofstream idf(prefix + ".grm.id");
  if (!idf) die("cannot write " + prefix + ".grm.id");
  for (const auto& id : g.ids) idf << id << '\t' << id << '\n';
  const uint64_t ntri = static_cast<uint64_t>(n) * (n + 1) / 2;
  std::vector<float> buf(ntri);
  size_t k = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j) buf[k++] = static_cast<float>(g.K(i, j));
  FILE* fp = std::fopen((prefix + ".grm.bin").c_str(), "wb");
  if (!fp) die("cannot write " + prefix + ".grm.bin");
  std::fwrite(buf.data(), sizeof(float), ntri, fp);
  std::fclose(fp);
  info("wrote GRM " + prefix + ".grm.id/.grm.bin");
}

Grm slice_grm(const Grm& g, const std::vector<std::string>& sample_ids) {
  auto m = index_map(g.ids);
  const int n = static_cast<int>(sample_ids.size());
  Grm o;
  o.ids = sample_ids;
  o.K.resize(n, n);
  std::vector<int> ix(n);
  for (int i = 0; i < n; ++i) {
    auto it = m.find(sample_ids[i]);
    if (it == m.end()) die("sample missing in GRM: " + sample_ids[i]);
    ix[i] = it->second;
  }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) o.K(i, j) = g.K(ix[i], ix[j]);
  return o;
}

namespace {

template <typename ForEach>
Grm compute_grm_from(ForEach&& for_each, const std::vector<std::string>& sample_ids) {
  // GRM audit vs GCTA --make-grm (src/GRM.cpp: GRM::calculate_GRM_blas,
  // PgenReader::CountHardDosage, Geno.cpp maf/min_maf/--geno):
  //  - standardization ((x - 2p) / sqrt(2 p (1 - p))) with p = ALT freq from
  //    non-missing alleles: identical to GCTA (mu = 2*af, sd = sqrt(2p(1-p))).
  //  - monomorphic rejection: eqtl p ∈ (1e-8, 1-1e-8); GCTA dev < 1e-50 — both
  //    drop effectively-monomorphic SNPs.
  //  - MAF filter: eqtl keeps maf_min <= min(af,1-af) (pass_maf); GCTA keeps
  //    maf >= min_maf with SMALL_EPSILON slack — equivalent in practice.
  //  - --max-miss default 0.8 (eqtl) vs --geno default 1.0 (GCTA): eqtl is
  //    strictly safer by default; set --max-miss 1.0 to match GCTA.
  //  - DIVISOR: eqtl divides every pair by the same m = #SNPs passing all
  //    marker filters. GCTA divides pair (i,j) by the count of SNPs jointly
  //    non-missing for (i,j) (src/GRM.cpp and sub_N in the BLAS write path).
  //    Because eqtl's parse_record mean-imputes missing dosage *before* the
  //    cross product, the off-diagonal numerator equals GCTA's (sum over
  //    jointly-observed SNPs); only the denominator differs. With complete
  //    genotypes (the common case for imputed/eQTL data) sub_N == m and the
  //    two are identical. With missing GT + --miss-hand impute this is a real
  //    bias for pairs where a sample is missing at >0 markers: eqtl
  //    under-estimates their relatedness. ALIGNMENT IS DEFERRED to a follow-up
  //    change that exposes a per-marker missing sample set (SnpRec miss mask)
  //    so denom_ij = m - miss_i - miss_j + miss_ij can be computed without
  //    breaking the LM/LMM/GLM callers that consume imputed dosage.
  const int n = static_cast<int>(sample_ids.size());
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  constexpr int kBatch = 256;
  Eigen::MatrixXd Z(n, kBatch); // ponytail: batch DGEMM instead of rank-1
  int col = 0, m = 0;
  for_each([&](const SnpRec& s) {
    Eigen::Map<const Eigen::VectorXd> g(s.dosage.data(), n);
    double mu = g.mean();
    double p = mu / 2.0;
    if (p < 1e-8 || p > 1.0 - 1e-8) return true;
    double sd = std::sqrt(2.0 * p * (1.0 - p));
    if (sd < 1e-12) return true;
    Z.col(col++) = (g.array() - mu) / sd;
    ++m;
    if (col == kBatch) {
      A.noalias() += Z * Z.transpose();
      col = 0;
    }
    return true;
  });
  if (col > 0) {
    A.noalias() += Z.leftCols(col) * Z.leftCols(col).transpose();
  }
  if (m == 0) die("no SNPs for GRM");
  A /= static_cast<double>(m);
  Grm grm;
  grm.ids = sample_ids;
  grm.K = std::move(A);
  info("computed GRM from " + std::to_string(m) + " SNPs");
  return grm;
}

}  // namespace

Grm compute_grm(VcfSession& vcf, const std::vector<std::string>& sample_ids, const MissPolicy& miss,
                double maf_min) {
  vcf.set_sample_order(sample_ids);
  return compute_grm_from(
      [&](auto&& fn) { vcf.for_each_snp(miss, maf_min, std::forward<decltype(fn)>(fn)); },
      sample_ids);
}

Grm compute_grm(PlinkBed& bed, const std::vector<std::string>& sample_ids, const MissPolicy& miss,
                double maf_min) {
  bed.set_sample_order(sample_ids);
  return compute_grm_from(
      [&](auto&& fn) { bed.for_each_snp(miss, maf_min, std::forward<decltype(fn)>(fn)); },
      sample_ids);
}

} // namespace eqtl
