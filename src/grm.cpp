#include "eqtl/grm.hpp"
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
  while (idf >> a) {
    // GCTA: FID IID or single ID
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

Grm compute_grm(VcfSession& vcf, const std::vector<std::string>& sample_ids,
                MissHand miss) {
  vcf.set_sample_order(sample_ids);
  const int n = static_cast<int>(sample_ids.size());
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, n);
  Eigen::VectorXd ones = Eigen::VectorXd::Ones(n);
  int m = 0;
  vcf.for_each_snp(miss, [&](const SnpRec& s) {
    Eigen::Map<const Eigen::VectorXd> g(s.dosage.data(), n);
    double mu = g.mean();
    double p = mu / 2.0;
    if (p < 1e-8 || p > 1.0 - 1e-8) return true;
    double sd = std::sqrt(2.0 * p * (1.0 - p));
    if (sd < 1e-12) return true;
    Eigen::VectorXd z = (g.array() - mu) / sd;
    A.noalias() += z * z.transpose();
    ++m;
    return true;
  });
  if (m == 0) die("no SNPs for GRM");
  A /= static_cast<double>(m);
  Grm grm;
  grm.ids = sample_ids;
  grm.K = std::move(A);
  info("computed GRM from " + std::to_string(m) + " SNPs");
  return grm;
}

} // namespace eqtl
