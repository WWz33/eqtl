#include "eqtl/pheno.hpp"
#include "eqtl/util.hpp"
#include <fstream>
#include <sstream>
#include <set>
#include <limits>
#include <cmath>
#include <cstdlib>

namespace eqtl {

static std::vector<std::string> split_line(const std::string& line) {
  // tab-preferred; fall back to whitespace
  if (line.find('\t') != std::string::npos) return split_char(line, '\t');
  return split_ws(line);
}

PhenoData load_pheno(const std::string& path) {
  std::ifstream in(path);
  if (!in) die("cannot open pheno: " + path);
  std::string line;
  if (!std::getline(in, line)) die("empty pheno: " + path);
  line = trim(line);
  if (!line.empty() && line[0] == '#') line = line.substr(1);
  auto header = split_line(line);
  if (header.size() < 2) die("pheno header needs sample + >=1 gene");
  // header[0] is sample column name (ignored); rest gene ids
  PhenoData P;
  for (size_t i = 1; i < header.size(); ++i) P.gene_ids.push_back(trim(header[i]));

  std::vector<std::vector<double>> rows;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    auto tok = split_line(line);
    if (tok.size() != header.size())
      die("pheno row column mismatch for sample " + (tok.empty() ? "?" : tok[0]));
    P.sample_ids.push_back(trim(tok[0]));
    std::vector<double> r(P.gene_ids.size());
    for (size_t j = 0; j < P.gene_ids.size(); ++j) {
      if (tok[j + 1] == "NA" || tok[j + 1] == "NaN" || tok[j + 1] == ".")
        r[j] = std::numeric_limits<double>::quiet_NaN();
      else r[j] = std::stod(tok[j + 1]);
    }
    rows.push_back(std::move(r));
  }
  if (P.sample_ids.empty()) die("no samples in pheno");
  // check dup ids
  {
    std::set<std::string> s(P.sample_ids.begin(), P.sample_ids.end());
    if (s.size() != P.sample_ids.size()) die("duplicate sample ids in pheno");
  }
  P.Y.resize(P.sample_ids.size(), P.gene_ids.size());
  for (size_t i = 0; i < rows.size(); ++i)
    for (size_t j = 0; j < P.gene_ids.size(); ++j) P.Y(i, j) = rows[i][j];
  info("pheno: " + std::to_string(P.sample_ids.size()) + " samples, " +
       std::to_string(P.gene_ids.size()) + " genes");
  return P;
}

void slice_pheno(PhenoData& p, const std::vector<std::string>& sample_ids) {
  auto m = index_map(p.sample_ids);
  Eigen::MatrixXd Y(sample_ids.size(), p.gene_ids.size());
  for (size_t i = 0; i < sample_ids.size(); ++i) {
    auto it = m.find(sample_ids[i]);
    if (it == m.end()) die("internal: sample missing in pheno slice");
    Y.row(i) = p.Y.row(it->second);
  }
  p.sample_ids = sample_ids;
  p.Y = std::move(Y);
}

CovData load_covar(const std::string& path, const std::vector<std::string>& sample_ids) {
  CovData C;
  C.sample_ids = sample_ids;
  if (path.empty()) {
    C.cov_names = {"intercept"};
    C.X = Eigen::MatrixXd::Ones(sample_ids.size(), 1);
    return C;
  }
  std::ifstream in(path);
  if (!in) die("cannot open covar: " + path);
  std::string line;
  std::vector<std::string> file_ids;
  std::vector<std::vector<double>> rows;
  std::vector<std::string> names;
  bool first = true;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    auto tok = split_line(line);
    if (tok.size() < 2) continue;
    // Header: col0 is sample id name (non-numeric) OR second token non-numeric.
    // Purely numeric second field ⇒ data row (even if col0 looks like text sample id).
    if (first) {
      first = false;
      auto is_num = [](const std::string& s) -> bool {
        if (s.empty()) return false;
        char* end = nullptr;
        std::strtod(s.c_str(), &end);
        return end && *end == '\0';
      };
      const bool c0_num = is_num(tok[0]);
      const bool c1_num = is_num(tok[1]);
      // header if second col not a number, or first col is "sample"/"id"-like non-number
      // and all remaining tokens are non-numeric names
      if (!c1_num || (!c0_num && !c1_num)) {
        names.assign(tok.begin() + 1, tok.end());
        continue;
      }
    }
    file_ids.push_back(tok[0]);
    std::vector<double> r;
    for (size_t j = 1; j < tok.size(); ++j) r.push_back(std::stod(tok[j]));
    rows.push_back(std::move(r));
  }
  if (rows.empty()) die("no covariates in " + path);
  int nc = static_cast<int>(rows[0].size());
  if (names.empty()) {
    for (int j = 0; j < nc; ++j) names.push_back("cov" + std::to_string(j + 1));
  }
  auto fmap = index_map(file_ids);
  Eigen::MatrixXd raw(sample_ids.size(), nc);
  for (size_t i = 0; i < sample_ids.size(); ++i) {
    auto it = fmap.find(sample_ids[i]);
    if (it == fmap.end()) die("sample " + sample_ids[i] + " missing in covar");
    if (static_cast<int>(rows[it->second].size()) != nc) die("covar row width mismatch");
    for (int j = 0; j < nc; ++j) raw(static_cast<int>(i), j) = rows[it->second][j];
  }
  bool has_int = false;
  for (int j = 0; j < nc; ++j) {
    bool ok = true;
    for (int i = 0; i < static_cast<int>(sample_ids.size()); ++i)
      if (std::fabs(raw(i, j) - 1.0) > 1e-12) {
        ok = false;
        break;
      }
    if (ok) {
      has_int = true;
      break;
    }
  }
  if (!has_int) {
    C.X.resize(sample_ids.size(), nc + 1);
    C.X.col(0).setOnes();
    C.X.rightCols(nc) = raw;
    C.cov_names.push_back("intercept");
    for (auto& n : names) C.cov_names.push_back(n);
  } else {
    C.X = raw;
    C.cov_names = names;
  }
  info("covar: " + std::to_string(C.X.cols()) + " columns (incl intercept if added)");
  return C;
}

} // namespace eqtl
