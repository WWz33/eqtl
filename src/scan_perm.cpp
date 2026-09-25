#include "eqtl/scan_trans.hpp"
#include <atomic>
#include <limits>

namespace eqtl {

// The trans stage-2 permutation is disabled. It compared the observed
// genome-wide min-p against a permutation null built on a top-K SNP set that
// was itself selected on the observed p-values: the observed side carries a
// selection bias the permutation side does not have, so p_emp pinned at
// 1/(B+1) for ~97% of genes on the null panel (measured; analysis and the
// replacement design in eqtl-gene-level-p-plan.md). Report NA rather than a
// number that reads as significant for every gene. cis permutations
// (scan_cis.hpp) are unaffected — their null covers the whole window with no
// p-based selection.
template <typename Job>
void stage2_perm_topk(const Options& opt, Model /*model*/, Job& job,
                      const LmmBasis* /*ext_basis*/) {
  job.summary.p_emp = std::numeric_limits<double>::quiet_NaN();
  job.summary.p_beta = std::numeric_limits<double>::quiet_NaN();
  static std::atomic<int> warned{0};
  if (opt.perm > 0 && !warned.exchange(1))
    warn("trans: stage-2 permutation disabled — its top-K null is biased "
         "(p_emp was pinned at 1/(B+1)); p_emp/p_beta are NA");
}

// Explicit instantiations
template void stage2_perm_topk<GeneLmJob>(const Options&, Model, GeneLmJob&, const LmmBasis*);
template void stage2_perm_topk<GeneLmmJob>(const Options&, Model, GeneLmmJob&, const LmmBasis*);

} // namespace eqtl
