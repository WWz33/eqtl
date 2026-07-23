/* eqtl — phenotype / covariates */
#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>

namespace eqtl {

struct PhenoData {
  std::vector<std::string> sample_ids; // analysis order
  std::vector<std::string> gene_ids;
  // samples x genes
  Eigen::MatrixXd Y;
};

struct CovData {
  std::vector<std::string> sample_ids;
  std::vector<std::string> cov_names;
  // samples x cov (includes intercept as col0 if added)
  Eigen::MatrixXd X;
};

// Load pheno: header row gene names; col0 = sample id; rest = values.
// sample_order preference = file order; will be intersected later.
PhenoData load_pheno(const std::string& path);

// Optional covar: col0 sample id. Adds intercept if none detected.
// Empty path -> intercept-only for given sample_ids.
CovData load_covar(const std::string& path, const std::vector<std::string>& sample_ids);

// Restrict pheno columns/rows to sample_ids order.
void slice_pheno(PhenoData& p, const std::vector<std::string>& sample_ids);

} // namespace eqtl
