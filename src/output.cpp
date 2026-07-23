#include "eqtl/output.hpp"
#include "eqtl/util.hpp"
#include <cstdio>
#include <cmath>

namespace eqtl {

void write_pairs_header(std::ostream& os, Model m) {
  os << "gene\tsnp\tchrom\tpos\tref\talt\tmaf\tbeta\tse\tstat\tp\tr2\tn\ttss_dist\tscope";
  if (m == Model::Glm) os << "\tphi\tglm_converged";
  if (m == Model::Glmm) os << "\tglmm_converged";
  os << '\n';
}

void write_pair_line(std::ostream& os, const AssocHit& h, Model m, const std::string& scope) {
  char buf[1024];
  int n;
  if (h.has_tss_dist) {
    n = std::snprintf(
        buf, sizeof(buf),
        "%s\t%s\t%s\t%lld\t%s\t%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%d\t%.10g\t%s",
        h.gene.c_str(), h.snp.c_str(), h.chrom.c_str(), static_cast<long long>(h.pos), h.ref.c_str(),
        h.alt.c_str(), h.maf, h.beta, h.se, h.stat, h.p, h.r2, h.n, h.tss_dist, scope.c_str());
  } else {
    n = std::snprintf(
        buf, sizeof(buf),
        "%s\t%s\t%s\t%lld\t%s\t%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%d\tNA\t%s",
        h.gene.c_str(), h.snp.c_str(), h.chrom.c_str(), static_cast<long long>(h.pos), h.ref.c_str(),
        h.alt.c_str(), h.maf, h.beta, h.se, h.stat, h.p, h.r2, h.n, scope.c_str());
  }
  if (n < 0) return;
  if (m == Model::Glm) {
    const int n2 = std::snprintf(buf + n, sizeof(buf) - static_cast<size_t>(n), "\t%.10g\t%d", h.phi,
                                 h.glm_converged ? 1 : 0);
    if (n2 > 0) n += n2;
  }
  if (m == Model::Glmm) {
    const int n2 =
        std::snprintf(buf + n, sizeof(buf) - static_cast<size_t>(n), "\t%d", h.glmm_converged ? 1 : 0);
    if (n2 > 0) n += n2;
  }
  if (n > 0 && static_cast<size_t>(n) < sizeof(buf) - 1) {
    buf[n++] = '\n';
    os.write(buf, n);
    return;
  }
  // long ids: stream fallback
  os << h.gene << '\t' << h.snp << '\t' << h.chrom << '\t' << h.pos << '\t' << h.ref << '\t'
     << h.alt << '\t' << h.maf << '\t' << h.beta << '\t' << h.se << '\t' << h.stat << '\t' << h.p
     << '\t' << h.r2 << '\t' << h.n << '\t';
  if (h.has_tss_dist) os << h.tss_dist;
  else os << "NA";
  os << '\t' << scope;
  if (m == Model::Glm) os << '\t' << h.phi << '\t' << (h.glm_converged ? 1 : 0);
  if (m == Model::Glmm) os << '\t' << (h.glmm_converged ? 1 : 0);
  os << '\n';
}

void write_top_header(std::ostream& os, Model m) { write_pairs_header(os, m); }

void write_region_header(std::ostream& os) {
  os << "gene\tchrom\ttss\tn_tested\tn_sig\tacat_p\tp_emp\tp_beta\tbeta_shape1\tbeta_shape2\n";
}

static void fmt_na_or(char* out, size_t n, double v) {
  if (!std::isfinite(v)) {
    std::snprintf(out, n, "NA");
  } else {
    std::snprintf(out, n, "%.10g", v);
  }
}

void write_region_line(std::ostream& os, const GeneSummary& g) {
  char p_emp[32], p_beta[32], buf[512];
  fmt_na_or(p_emp, sizeof(p_emp), g.p_emp);
  fmt_na_or(p_beta, sizeof(p_beta), g.p_beta);
  const int n = std::snprintf(
      buf, sizeof(buf), "%s\t%s\t%lld\t%d\t%d\t%.10g\t%s\t%s\t%.10g\t%.10g\n", g.gene.c_str(),
      g.chrom.c_str(), static_cast<long long>(g.tss), g.n_tested, g.n_sig, g.acat_p, p_emp, p_beta,
      g.beta_shape1, g.beta_shape2);
  if (n > 0) os.write(buf, n);
}

}  // namespace eqtl
