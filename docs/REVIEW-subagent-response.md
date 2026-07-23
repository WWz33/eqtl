# Subagent review response (deleg_1bf96b74) — updated

## Closed vs open

| Review item | Status |
|-------------|--------|
| LMM eig every gene | **Fixed** (`make_lmm_basis` once) |
| cis full-VCF per gene | **Fixed** (cis-only: TBI/CSI region query; all/trans: load once) |
| residualise+score | never shipped |
| GT only / four models / CLI shape | OK |
| LMM logdet via determinant | **Fixed** (LDLT `vectorD` log-sum) |
| `--thread` unused | **Fixed** (OpenMP over SNPs within gene) |
| seed=-1 fixed to 1 | **Fixed** (`random_device`) |
| help missing `gw` | **Fixed** |
| `for_each_snp_region` stub | **Fixed** (CSI/`bcf_itr` or TBI/`tbx_itr`; fallback sequential) |
| chrom `chr1` vs `1` | **Fixed** (`chrom_equal` / `chrom_key`) |
| miss filter too strict only | **Improved** (`--max-miss FLOAT` [0]) |
| lm numeric gold | **Added** `scripts/gold_lm.py` (numpy OLS; smoke row ≤1e-4 beta/se) |
| `--perm` default cost | **Open** (skeleton works; default 10000 still heavy) |
| pheno per-sample NA | **Open** |
| GEMMA lmm gold | **Open** |
| vendored htslib default | **Open** (`USE_SYSTEM_HTS=1`) |
| gene-level OpenMP + thread-local writes | **Open** (SNP-parallel only) |

## Verify (ad-hoc)

`AD-HOC VERIFY (review-continue): ALL PASS` — version, help, cis region path, max-miss, grm+lmm eig once, mode all load, gold_lm.
