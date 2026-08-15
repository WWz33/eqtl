/* eqtl — shared scan types and per-SNP test dispatch */
#pragma once
#include "eqtl/options.hpp"
#include "eqtl/pheno.hpp"
#include "eqtl/annot.hpp"
#include "eqtl/grm.hpp"
#include "eqtl/models.hpp"
#include "eqtl/output.hpp"
#include "eqtl/stats_extra.hpp"
#include "eqtl/util.hpp"
#include "eqtl/vcf_session.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <cmath>
#include <utility>

namespace eqtl {

// ---------------------------------------------------------------------------
// Output file set for one (model, scope) combination
// ---------------------------------------------------------------------------
struct ScopeOut {
  std::ofstream pairs;
  std::ofstream top;
  std::ofstream region;
  std::string tag;
  std::vector<char> pairs_buf;
  std::vector<char> top_buf;
  std::vector<char> region_buf;
};

// ---------------------------------------------------------------------------
// Per-gene prepared data (subset to complete observations)
// ---------------------------------------------------------------------------
struct GeneReady {
  std::vector<int> keep;
  Eigen::VectorXd y;
  Eigen::MatrixXd X;
  Eigen::MatrixXd K;
  LmmBasis basis;
  bool has_basis = false;
  // Optional shared refs: when set, prep_null uses these instead of copying K/basis.
  const Eigen::MatrixXd* K_ref = nullptr;
  const LmmBasis* basis_ref = nullptr;
};

// ---------------------------------------------------------------------------
// Small inline helpers
// ---------------------------------------------------------------------------
inline bool needs_grm(Model m) { return m == Model::Lmm || m == Model::Glmm; }
inline bool needs_counts(Model m) { return m == Model::Glm || m == Model::Glmm; }

// cis window anchored on the gene body [start-window, end+window], matching
// QTLtools / GTEx / eQTLcatalogue; a superset of TSS±window. window is the
// flanking distance in bp (phenotype_start/phenotype_end are 1-based inclusive).
template <typename T>
inline std::pair<int64_t, int64_t> cis_window_bounds(const GeneLoc& loc, T window) {
  int64_t span_start = loc.start - static_cast<int64_t>(window);
  if (span_start < 1) span_start = 1;
  return {span_start, loc.end + static_cast<int64_t>(window)};
}

inline bool in_cis_window(const SnpRec& s, const GeneLoc& loc, int window) {
  if (!chrom_equal(s.chrom, loc.chrom)) return false;
  const auto [cstart, cend] = cis_window_bounds(loc, window);
  return s.pos >= cstart && s.pos <= cend;
}

// ---------------------------------------------------------------------------
// Gene preparation
// ---------------------------------------------------------------------------
bool build_gene_ready(const Eigen::VectorXd& y_full, const Eigen::MatrixXd& X_full,
                      const Eigen::MatrixXd* K_full, bool need_k, bool need_lmm_basis,
                      bool fast_sparse, GeneReady& out);

// Shared-basis fast path skips build_gene_ready's per-sample filtering, so it
// requires y and every covariate column to be all-finite. Non-finite entries
// (Inf from overflow, etc.) would otherwise reach Q^T*X and silently fail LDLT.
inline bool all_finite(const Eigen::VectorXd& y, const Eigen::MatrixXd& X) {
  for (Eigen::Index i = 0; i < y.size(); ++i)
    if (!std::isfinite(y(i))) return false;
  for (Eigen::Index j = 0; j < X.cols(); ++j)
    for (Eigen::Index i = 0; i < X.rows(); ++i)
      if (!std::isfinite(X(i, j))) return false;
  return true;
}

Eigen::VectorXd subset_dosage(const std::vector<double>& full, const std::vector<int>& keep);

// ---------------------------------------------------------------------------
// Per-SNP test dispatch (model-generic)
// ---------------------------------------------------------------------------
AssocHit run_test(Model model, bool fast, const GeneReady& gr, const Eigen::VectorXd& g,
                  GenePrepLm* lm_cache, GenePrepLmm* lmm_cache, GenePrepGlm* glm_cache,
                  GenePrepGlmm* glmm_cache, bool have_cache);

void prep_null(Model model, bool fast, const GeneReady& gr, GenePrepLm* lm_cache,
               GenePrepLmm* lmm_cache, GenePrepGlm* glm_cache, GenePrepGlmm* glmm_cache);

double subset_maf_or_nan(const Eigen::VectorXd& g, double* maf_out);

void fill_snp_id(AssocHit& h, const SnpRec& snp);

AssocHit test_one(Model model, bool fast, const GeneReady& gr, const SnpRec& snp,
                  const std::string& gene, const GeneLoc* loc, GenePrepLm* lm_c, GenePrepLmm* lmm_c,
                  GenePrepGlm* glm_c, GenePrepGlmm* glmm_c, bool have_cache,
                  Eigen::VectorXd& g_buf);

double test_one_p(Model model, bool fast, const GeneReady& gr, const Eigen::VectorXd& g,
                  GenePrepLm* lm_c, GenePrepLmm* lmm_c, GenePrepGlm* glm_c, GenePrepGlmm* glmm_c);

// ---------------------------------------------------------------------------
// Multiple testing
// ---------------------------------------------------------------------------
void bh_adjust(std::vector<GeneSummary>& gs);

} // namespace eqtl
