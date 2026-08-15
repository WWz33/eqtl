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
  // emit per-pair denominator count (.grm.N.bin, float32 little-incl-diag) so
  // the produced files are interoperable with GCTA's downstream (e.g. --grim-
  // reml reads .grm.N.bin). Computed only by compute_grm; absent N_packed ⇒
  // emit map-wide denominator as a uniform m-like alternative — but current
  // callers always populate N_packed. If empty we skip the file silently,
  // matching the historical output for unknown-denominator GRMs (load path).
  if (!g.N_packed.empty()) {
    if (g.N_packed.size() != ntri) die("GRM N_packed size mismatch on write");
    for (size_t i = 0; i < ntri; ++i) buf[i] = static_cast<float>(g.N_packed[i]);
    FILE* nfp = std::fopen((prefix + ".grm.N.bin").c_str(), "wb");
    if (!nfp) die("cannot write " + prefix + ".grm.N.bin");
    std::fwrite(buf.data(), sizeof(float), ntri, nfp);
    std::fclose(nfp);
    info("wrote GRM " + prefix + ".grm.id/.grm.bin/.grm.N.bin");
  } else {
    info("wrote GRM " + prefix + ".grm.id/.grm.bin");
  }
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
  // reorder per-pair denominator (n*(n+1)/2 packed lower-incl-diag) the same way
  // as K so write_grm_gcta still emits matching .grm.N.bin; only if present.
  if (!g.N_packed.empty()) {
    o.N_packed.resize(static_cast<size_t>(n) * (n + 1) / 2);
    auto get_N = [&](int I, int J) {
      if (I < J) std::swap(I, J);
      return g.N_packed[static_cast<size_t>(I) * (I + 1) / 2 + J];
    };
    size_t k = 0;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j <= i; ++j)
        o.N_packed[k++] = get_N(ix[i], ix[j]);
  }
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
  //  - MISSING-IMPUTATION vs BLAS X: GCTA getGenoDouble_bed sets standardized
  //    X = (geno - center_value) * rdev with center_value = mu; missing calls
  //    take the "na" lookup entry = (mu - mu) * rdev = 0 (Geno.cpp ~951-974).
  //    eqtl mean-imputes missing dosage with the non-missing sample mean,
  //    giving Z_missing = (mu - mu)/sd = 0 — identical numerator contribution.
  //  - DIVISOR: GCTA divides pair (i,j) by sub_N(i,j) = numValidMarkers −
  //    sub_miss[i] − sub_miss[j] + po_N(i,j), where sub_miss[k] is the total
  //    marker-missing count for sample k and po_N(i,j) is the jointly-missing
  //    count (src/GRM.cpp:1442-1448). Because mean-imputation zeroes out the
  //    missing-sample cross-product, eqtl's accumulated numerator A[i,j] =
  //    Σ_{markers jointly-observed for (i,j)} Z_i Z_j already equals GCTA's; we
  //    only had the divisor wrong (m, not sub_N). We now track per-sample
  //    sub_miss and lower-triangular joint_miss from each SnpRec's miss_mask
  //    (captured before imputation by vcf_session/plink_bed) and divide each
  //    pair by sub_N. Under --miss-hand filter (no missing passes), sub_miss
  //    and joint_miss stay 0 and sub_N == m — identical to the old behavior
  //    and to GCTA's no-missing case.
  const int n = static_cast<int>(sample_ids.size());
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  // per-sample marker-missing count and pair-jointly-missing count (lower tri,
  // diagonal implicit; po_N(i,i) == sub_miss[i] in GCTA's formula).
  Eigen::VectorXi sub_miss = Eigen::VectorXi::Zero(n);
  Eigen::MatrixXi joint_miss = Eigen::MatrixXi::Zero(n, n);
  constexpr int kBatch = 256;
  Eigen::MatrixXd Z(n, kBatch); // ponytail: batch DGEMM instead of rank-1
  int col = 0, m = 0;
  std::vector<int> miss_idx;   // per-SNP missing sample indices (reused buffer)
  miss_idx.reserve(static_cast<size_t>(n));
  for_each([&](const SnpRec& s) {
    Eigen::Map<const Eigen::VectorXd> g(s.dosage.data(), n);
    double mu = g.mean();
    double p = mu / 2.0;
    if (p < 1e-8 || p > 1.0 - 1e-8) return true;
    double sd = std::sqrt(2.0 * p * (1.0 - p));
    if (sd < 1e-12) return true;
    Z.col(col) = (g.array() - mu) / sd;  // missing → 0 (mean-imputed before here)
    ++col;
    ++m;
    // Cumulative per-pair non-missing denominator = m − union(missing_i, missing_j).
    // Track per-sample sub_miss[k] and lower-tri joint_miss(i,j) for the
    // missing samples of this SNP. Empty miss_mask (most SNPs) is O(1); only
    // SNPs with missingness pay the O(k²) inner loop, k = #missing.
    if (!s.miss_mask.empty()) {
      miss_idx.clear();
      for (int i = 0; i < n; ++i) {
        if (s.miss_mask[static_cast<size_t>(i)] != 0) {
          ++sub_miss[i];
          // j ∈ miss_idx always precedes i (iteration order) → store lower-tri
          for (const int j : miss_idx) joint_miss(i, j) += 1;
          miss_idx.push_back(i);
        }
      }
    }
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
  // Convert from "mean-imputed Z outer product summed over markers" to GCTA's
  // per-pair jointly-non-missing relatedness: A[i,j] / sub_N(i,j), sub_N = m −
  // miss_i − miss_j + joint_miss(i,j); diagonal sub_N(i,i) = m − sub_miss[i]
  // (= m − 2*sub_miss[i] + sub_miss[i]). Use the symmetric lower-triangular
  // joint_miss for both (i,j) and (j,i). With zero missingness (the filter
  // path or complete-genotype VCFs/beds), this collapses to dividing by m and
  // is binary-identical to the previous behavior.
  for (int i = 0; i < n; ++i) {
    const int mi = sub_miss[i];
    for (int j = 0; j < i; ++j) {
      // joint_miss stored strictly lower-tri (max_idx, min_idx)
      const int jm = joint_miss(i, j);
      const int sub_N = m - mi - sub_miss[j] + jm;
      double v = sub_N > 0 ? A(i, j) / static_cast<double>(sub_N) : 0.0;
      A(i, j) = v;
      A(j, i) = v;
    }
    const int sub_N_diag = m - mi;  // po_N(i,i) == sub_miss[i]
    A(i, i) = sub_N_diag > 0 ? A(i, i) / static_cast<double>(sub_N_diag) : 0.0;
  }
  Grm grm;
  grm.ids = sample_ids;
  grm.K = std::move(A);
  // emit per-pair jointly-non-missing marker count (GCTA .grm.N.bin): packed
  // n*(n+1)/2, row-major lower-incl-diag (idx = i*(i+1)/2 + j).
  grm.N_packed.resize(static_cast<size_t>(n) * (n + 1) / 2);
  size_t k = 0;
  for (int i = 0; i < n; ++i) {
    const int mi = sub_miss[i];
    for (int j = 0; j < i; ++j) {
      const int sub_N = m - mi - sub_miss[j] + joint_miss(i, j);
      grm.N_packed[k++] = sub_N;
    }
    grm.N_packed[k++] = m - mi;  // po_N(i,i) == sub_miss[i]
  }
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
