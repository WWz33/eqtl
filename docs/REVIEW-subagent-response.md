# Subagent review response (deleg_1bf96b74)

Review timestamp was ~04:08; several items were already fixed in `be4719c` before/around review.

## Already fixed (before or independent of this batch)

| Review critical | Status in tree now |
|-----------------|-------------------|
| LMM eig every gene | **Fixed**: `make_lmm_basis` once in `run_eqtl` |
| cis full-VCF per gene | **Fixed**: load `all_snps` once; cis/trans filter in memory |
| residualise+score product | **OK** (never shipped) |
| GT only, models, CLI shape | **OK** |

## Fixed in this follow-up

1. LMM REML `log|X'D^{-1}X|` via LDLT `vectorD()` log-sum (not `determinant()`)
2. `--thread`: OpenMP over **SNPs within gene** (prep once; write serial)
3. `--seed < 0` → `std::random_device` (not fixed seed 1)
4. help/README: mode documents `gw`

## Still open (v1.1)

1. Real TBI/CSI region iterators (currently sequential load-all; fine for smoke, not for huge VCF memory)
2. `--perm` cost (default 10000 full refits)
3. miss-rate threshold; pheno per-sample NA
4. chrom name normalize (`chr1` vs `1`)
5. Numeric gold vs MatrixEQTL/GEMMA
6. Vendored htslib default
7. Gene-level OpenMP with thread-local buffers (current: SNP-parallel only)

## Note on review item “for_each_snp_region unused”

Still true as API stub; product path no longer calls it for cis (uses in-memory `all_snps`). Index-backed streaming remains v1.1 work.
