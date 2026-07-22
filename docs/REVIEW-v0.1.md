# eqtl v0.1 — first-pass review (local)

Date: 2026-07-23  
Tree: `/home/ww/eqtl/eqtl`  
Binary: `eqtl 0.1.0`

## Scope

Product surface + numeric engine vs LOCK (grill session). Not a full GEMMA/MatrixEQTL gold campaign.

## Critical (fixed in this pass)

1. **LMM re-eigensolved GRM every gene** (`prep_lmm` called SelfAdjointEigenSolver on K each prep).  
   - Fix: `make_lmm_basis(K)` once in `run_eqtl`; `prep_lmm(y,X,basis,fast)`.

2. **cis re-scanned whole VCF per gene** (`collect_cis` → sequential full-file filter).  
   - Fix: load SNPs once with `load_all`; cis/trans window in memory.

## Critical / high (open)

1. **No automated numeric gold** vs MatrixEQTL (lm) or GEMMA `-lmm` (≤1e-4). Smoke only checks non-empty outputs.
2. **`--perm` default 10000 is extremely expensive** (full model refit per SNP × B). Usable only with `--perm 0` for scans; need FastQTL-style residual reuse or subset.
3. **`--miss-hand filter` default drops any SNP with any missing GT** among analysis samples — often near-empty real panels; document and consider fraction threshold later.
4. **glmm / NB IRLS** are self-written; not yet bit/float-aligned to a reference implementation.

## Medium

1. VCF index: `bcf_index_load` may fail on .vcf.gz while `tbx_index_load` works; stderr noise from htslib; region TBI path still simplified to sequential filter in `for_each_snp_region`.
2. Multi-allelic SNPs skipped; no DS path (LOCK OK).
3. `load_grm_gcta` FID/IID parsing is whitespace-fragile for odd ID files.
4. OpenMP flag set but gene loop still serial (no `#pragma omp` yet).
5. `third_party/eigen` vendored full tree; htslib still system default (`USE_SYSTEM_HTS=1`).

## Low / style

1. Help/README match mainstream-bio-cli-style (behavior-only) — OK for v0.1.
2. Dead `collect_cis` removed after in-memory cis.
3. Temporary smoke outputs under `data/*_*.tsv` not gitignored patterns for all prefixes — only `data/smoke_out*` in ignore.

## LOCK compliance (summary)

| Item | Status |
|------|--------|
| Single CLI `eqtl`, `--make-grm` | OK |
| Models lm/glm/lmm/glmm | OK (paths present) |
| Not quasar residualise+score | OK |
| GT only | OK |
| GFF TSS, modes, outputs TSV | OK |
| emp p + ACAT | OK (perm costly) |
| `--fast` | Wired in prep/test |
| Numeric gold | NOT done |

## Verify (ad-hoc)

See session ad-hoc script: make, help, make-grm, lm cis, lmm cis — passed before this optimization pass. Re-run after rebuild.

## v1.1 suggested order

1. Rebuild + smoke after LMM basis / SNP load fix  
2. Tiny lm gold vs numpy/MatrixEQTL  
3. Tiny lmm gold vs GEMMA  
4. Soften miss filter or document impute as recommended  
5. Speed: OpenMP over genes; faster perm  
6. Vendor htslib submodule for clone-build
