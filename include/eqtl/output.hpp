/* eqtl — output writers */
#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include "eqtl/models.hpp"
#include "eqtl/options.hpp"

namespace eqtl {

struct GeneSummary {
  std::string gene;
  std::string chrom;
  int64_t tss = 0;
  double acat_p = 1;
  double p_emp = std::numeric_limits<double>::quiet_NaN();
  double p_beta = std::numeric_limits<double>::quiet_NaN();
  double beta_shape1 = 0;
  double beta_shape2 = 0;
  int n_tested = 0;
  int n_sig = 0;
};

void write_pairs_header(std::ostream& os, Model m);
void write_pair_line(std::ostream& os, const AssocHit& h, Model m, const std::string& scope);
void write_top_header(std::ostream& os, Model m);
void write_region_header(std::ostream& os);
void write_region_line(std::ostream& os, const GeneSummary& g);

} // namespace eqtl
