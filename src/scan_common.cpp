#include "eqtl/scan_common.hpp"
#include <algorithm>
#include <numeric>
#include <cstring>

namespace eqtl {

// ---------------------------------------------------------------------------
// Gene preparation
// ---------------------------------------------------------------------------
bool build_gene_ready(const Eigen::VectorXd& y_full, const Eigen::MatrixXd& X_full,
                      const Eigen::MatrixXd* K_full, bool need_k, bool need_lmm_basis,
                      bool fast_sparse, GeneReady& out) {
  const int n_full = static_cast<int>(y_full.size());
  out.keep.clear();
  out.keep.reserve(static_cast<size_t>(n_full));
  for (int i = 0; i < n_full; ++i) {
    if (!std::isfinite(y_full(i))) continue;
    bool ok = true;
    for (int j = 0; j < X_full.cols(); ++j) {
      if (!std::isfinite(X_full(i, j))) { ok = false; break; }
    }
    if (ok) out.keep.push_back(i);
  }
  const int nk = static_cast<int>(out.keep.size());
  if (nk < 3) return false;

  out.y.resize(nk);
  out.X.resize(nk, X_full.cols());
  for (int r = 0; r < nk; ++r) {
    const int i = out.keep[static_cast<size_t>(r)];
    out.y(r) = y_full(i);
    out.X.row(r) = X_full.row(i);
  }
  if (out.y.array().abs().maxCoeff() < 1e-15) return false;

  out.K.resize(0, 0);
  out.K_ref = nullptr;
  out.basis_ref = nullptr;
  out.has_basis = false;
  if (need_k) {
    if (!K_full) return false;
    if (nk == n_full) {
      // keep is the full identity set → reference the shared GRM, no O(n²) copy
      out.K_ref = K_full;
    } else {
      out.K.resize(nk, nk);
      for (int a = 0; a < nk; ++a)
        for (int b = 0; b < nk; ++b)
          out.K(a, b) = (*K_full)(out.keep[static_cast<size_t>(a)], out.keep[static_cast<size_t>(b)]);
    }
    if (need_lmm_basis) {
      const Eigen::MatrixXd& K_src = out.K_ref ? *out.K_ref : out.K;
      if (fast_sparse) {
        Eigen::MatrixXd K_use = K_src;
        sparsify_grm(K_use, 1e-4);
        out.basis = make_lmm_basis(K_use);
      } else {
        out.basis = make_lmm_basis(K_src);
      }
      out.has_basis = true;
    }
  }
  return true;
}

Eigen::VectorXd subset_dosage(const std::vector<double>& full, const std::vector<int>& keep) {
  Eigen::VectorXd g(static_cast<int>(keep.size()));
  for (size_t r = 0; r < keep.size(); ++r) {
    const int i = keep[r];
    if (i < 0 || static_cast<size_t>(i) >= full.size())
      g(static_cast<int>(r)) = std::numeric_limits<double>::quiet_NaN();
    else
      g(static_cast<int>(r)) = full[static_cast<size_t>(i)];
  }
  return g;
}

// ---------------------------------------------------------------------------
// Per-SNP test dispatch
// ---------------------------------------------------------------------------
AssocHit run_test(Model model, bool fast, const GeneReady& gr, const Eigen::VectorXd& g,
                  GenePrepLm* lm_cache, GenePrepLmm* lmm_cache, GenePrepGlm* glm_cache,
                  GenePrepGlmm* glmm_cache, bool have_cache) {
  switch (model) {
    case Model::Lm: {
      if (!have_cache) *lm_cache = prep_lm(gr.y, gr.X);
      return test_lm(*lm_cache, g);
    }
    case Model::Lmm: {
      if (!have_cache) {
        if (gr.has_basis) *lmm_cache = prep_lmm(gr.y, gr.X, gr.basis, fast);
        else *lmm_cache = prep_lmm(gr.y, gr.X, gr.K_ref ? *gr.K_ref : gr.K, fast);
      }
      return test_lmm(*lmm_cache, g);
    }
    case Model::Glm: {
      if (!have_cache) *glm_cache = prep_glm_nb(gr.y, gr.X, fast);
      return test_glm_nb(*glm_cache, g);
    }
    case Model::Glmm: {
      if (!have_cache) *glmm_cache = prep_glmm_pois(gr.y, gr.X, gr.K_ref ? *gr.K_ref : gr.K, fast);
      return test_glmm_pois(*glmm_cache, g);
    }
  }
  return {};
}

void prep_null(Model model, bool fast, const GeneReady& gr, GenePrepLm* lm_cache,
               GenePrepLmm* lmm_cache, GenePrepGlm* glm_cache, GenePrepGlmm* glmm_cache) {
  switch (model) {
    case Model::Lm:
      *lm_cache = prep_lm(gr.y, gr.X);
      break;
    case Model::Lmm: {
      const LmmBasis& b = gr.basis_ref ? *gr.basis_ref : gr.basis;
      const Eigen::MatrixXd& k = gr.K_ref ? *gr.K_ref : gr.K;
      if (gr.has_basis) *lmm_cache = prep_lmm(gr.y, gr.X, b, fast);
      else *lmm_cache = prep_lmm(gr.y, gr.X, k, fast);
      break;
    }
    case Model::Glm:
      *glm_cache = prep_glm_nb(gr.y, gr.X, fast);
      break;
    case Model::Glmm: {
      const Eigen::MatrixXd& k = gr.K_ref ? *gr.K_ref : gr.K;
      *glmm_cache = prep_glmm_pois(gr.y, gr.X, k, fast);
      break;
    }
  }
}

double subset_maf_or_nan(const Eigen::VectorXd& g, double* maf_out) {
  const int n = static_cast<int>(g.size());
  if (n <= 0) return std::numeric_limits<double>::quiet_NaN();
  double sum = 0.0, sum2 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(g(i))) return std::numeric_limits<double>::quiet_NaN();
    sum += g(i);
    sum2 += g(i) * g(i);
  }
  double af = (sum / static_cast<double>(n)) / 2.0;
  if (af < 0) af = 0;
  if (af > 1) af = 1;
  if (maf_out) *maf_out = af;
  const double maf = af > 0.5 ? 1.0 - af : af;
  if (maf < 1e-12) return std::numeric_limits<double>::quiet_NaN();
  const double var = sum2 / n - (sum / n) * (sum / n);
  if (var < 1e-12) return std::numeric_limits<double>::quiet_NaN();
  return af;
}

void fill_snp_id(AssocHit& h, const SnpRec& snp) {
  if (!snp.id.empty()) h.snp = snp.id;
  else h.snp = snp.chrom + ":" + std::to_string(snp.pos) + ":" + snp.ref + ":" + snp.alt;
}

AssocHit test_one(Model model, bool fast, const GeneReady& gr, const SnpRec& snp,
                  const std::string& gene, const GeneLoc* loc, GenePrepLm* lm_c, GenePrepLmm* lmm_c,
                  GenePrepGlm* glm_c, GenePrepGlmm* glmm_c, bool have_cache,
                  Eigen::VectorXd& g_buf) {
  if (static_cast<int>(snp.dosage.size()) == static_cast<int>(gr.keep.size())) {
    if (g_buf.size() != static_cast<int>(snp.dosage.size()))
      g_buf.resize(static_cast<int>(snp.dosage.size()));
    std::memcpy(g_buf.data(), snp.dosage.data(), sizeof(double) * snp.dosage.size());
  } else {
    g_buf = subset_dosage(snp.dosage, gr.keep);
  }
  double maf_sub = snp.maf;
  if (!std::isfinite(subset_maf_or_nan(g_buf, &maf_sub))) {
    AssocHit h;
    h.p = std::numeric_limits<double>::quiet_NaN();
    return h;
  }
  AssocHit h = run_test(model, fast, gr, g_buf, lm_c, lmm_c, glm_c, glmm_c, have_cache);
  h.maf = maf_sub;
  h.n = static_cast<int>(gr.keep.size());
  h.gene = gene;
  h.chrom = snp.chrom;
  h.pos = snp.pos;
  h.ref = snp.ref;
  h.alt = snp.alt;
  if (loc && loc->ok) {
    h.has_tss_dist = true;
    h.tss_dist = static_cast<double>(snp.pos - loc->tss);
  }
  return h;
}

double test_one_p(Model model, bool fast, const GeneReady& gr, const Eigen::VectorXd& g,
                  GenePrepLm* lm_c, GenePrepLmm* lmm_c, GenePrepGlm* glm_c, GenePrepGlmm* glmm_c) {
  double maf_sub;
  if (!std::isfinite(subset_maf_or_nan(g, &maf_sub)))
    return std::numeric_limits<double>::quiet_NaN();
  return run_test(model, fast, gr, g, lm_c, lmm_c, glm_c, glmm_c, true).p;
}

// ---------------------------------------------------------------------------
// Multiple testing
// ---------------------------------------------------------------------------
void bh_adjust(std::vector<GeneSummary>& gs) {
  const int m = static_cast<int>(gs.size());
  if (m == 0) return;
  std::vector<int> idx(static_cast<size_t>(m));
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
            [&](int a, int b) { return gs[static_cast<size_t>(a)].acat_p < gs[static_cast<size_t>(b)].acat_p; });
  double prev = 1.0;
  for (int rank = m; rank >= 1; --rank) {
    auto& g = gs[static_cast<size_t>(idx[static_cast<size_t>(rank - 1)])];
    double q = g.acat_p * static_cast<double>(m) / static_cast<double>(rank);
    if (q > prev) q = prev;
    if (q > 1.0) q = 1.0;
    g.q_bh = q;
    prev = q;
  }
}

} // namespace eqtl
