/* eqtl — GRM */
#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "eqtl/vcf_session.hpp"
#include "eqtl/options.hpp"

namespace eqtl {

struct Grm {
  std::vector<std::string> ids;
  Eigen::MatrixXd K; // dense n x n
};

// Read GCTA .grm.id + .grm.bin (float32 lower incl diagonal).
Grm load_grm_gcta(const std::string& prefix);

// Write GCTA format.
void write_grm_gcta(const std::string& prefix, const Grm& g);

// Reorder/extract by sample_ids (error if missing).
Grm slice_grm(const Grm& g, const std::vector<std::string>& sample_ids);

// Compute GRM from VCF GT dosages on sample_order.
Grm compute_grm(VcfSession& vcf, const std::vector<std::string>& sample_ids, const MissPolicy& miss);

} // namespace eqtl
