#include "eqtl/output.hpp"
#include "eqtl/util.hpp"
#include <iomanip>
#include <sstream>
#include <cmath>

namespace eqtl {

void write_pairs_header(std::ostream& os, Model m) {
  os << "gene\tsnp\tchrom\tpos\tref\talt\tmaf\tbeta\tse\tstat\tp\tr2\tn\ttss_dist\tscope";
  if (m == Model::Glm) {
    os << "\tphi\tglm_converged";
  }
  if (m == Model::Glmm) {
    os << "\tglmm_converged";
  }
  os << '\n';
}

void write_pair_line(std::ostream& os, const AssocHit& h, Model m, const std::string& scope) {
  os << h.gene << '\t' << h.snp << '\t' << h.chrom << '\t' << h.pos << '\t' << h.ref << '\t'
     << h.alt << '\t' << h.maf << '\t' << h.beta << '\t' << h.se << '\t' << h.stat << '\t' << h.p
     << '\t' << h.r2 << '\t' << h.n << '\t';
  if (h.has_tss_dist) {
    os << h.tss_dist;
  } else {
    os << "NA";
  }
  os << '\t' << scope;
  if (m == Model::Glm) {
    os << '\t' << h.phi << '\t' << (h.glm_converged ? 1 : 0);
  }
  if (m == Model::Glmm) {
    os << '\t' << (h.glmm_converged ? 1 : 0);
  }
  os << '\n';
}

void write_top_header(std::ostream& os, Model m) {
  write_pairs_header(os, m);
}

void write_region_header(std::ostream& os) {
  os << "gene\tchrom\ttss\tn_tested\tn_sig\tacat_p\tp_emp\tp_beta\tbeta_shape1\tbeta_shape2\n";
}

void write_region_line(std::ostream& os, const GeneSummary& g) {
  auto na_or = [](double v) -> std::string {
    if (!std::isfinite(v)) return "NA";
    std::ostringstream ss;
    ss << v;
    return ss.str();
  };
  os << g.gene << '\t' << g.chrom << '\t' << g.tss << '\t' << g.n_tested << '\t' << g.n_sig << '\t'
     << g.acat_p << '\t' << na_or(g.p_emp) << '\t' << na_or(g.p_beta) << '\t' << g.beta_shape1
     << '\t' << g.beta_shape2 << '\n';
}

}  // namespace eqtl
