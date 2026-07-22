# eqtl

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

Map cis/trans expression QTLs from VCF genotypes and a sample×gene phenotype matrix.

## Getting Started

```bash
git clone https://github.com/WWz33/eqtl.git
cd eqtl && make -j

./scripts/make_smoke.sh
./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff \
  -c data/smoke.covar.tsv --model lm --mode cis --perm 0 --miss-hand impute -o data/out
```

## Usage

```text
eqtl [options]
```

| Flag | Default | Effect |
|------|---------|--------|
| `-v, --vcf` | required | VCF/BCF genotypes (GT→0/1/2) |
| `-e, --pheno` | required* | phenotype matrix (col1=sample) |
| `-g, --gff` | — | GFF3 for TSS; omit → genome-wide pairs |
| `-c, --covar` | — | covariates (optional) |
| `-k, --grm` | — | GCTA GRM prefix |
| `--make-grm` | off | write GRM from VCF and exit |
| `-m, --mode` | all | cis\|trans\|all |
| `--model` | lmm | lm\|glm\|lmm\|glmm (comma list) |
| `-w, --window` | 1000000 | cis window (bp) around TSS |
| `--pval-cis` | 1e-5 | cis write p threshold |
| `--pval-trans` | 1e-5 | trans/gw write p threshold |
| `--miss-hand` | filter | filter\|impute missing GT |
| `--fast` | off | per-gene shared nuisance |
| `--perm` | 10000 | gene-level permutations (0=off) |
| `--seed` | — | permutation seed |
| `--disable-beta-approx` | off | p_emp only |
| `-o, --out` | eqtl_out | output prefix |
| `-t, --thread` | 1 | threads |

\* not required with `--make-grm`

## Input / Output

| Kind | Notes |
|------|--------|
| Pheno | header gene names; col1 sample; remaining values |
| Covar | col1 sample id; intercept added if absent |
| GFF | gene features; TSS = start (+) or end (−) |
| GRM | `prefix.grm.id` + `prefix.grm.bin` (GCTA) |
| Pairs | `{out}.{model}.{cis\|trans\|gw}.pairs.tsv` |
| Top | `{out}.{model}.{scope}.top.tsv` |
| Region | `{out}.{model}.{scope}.region.tsv` (ACAT, p_emp, p_beta) |

## Citation

Software under development. Cite methods as appropriate for your analysis design.

## License

MIT
