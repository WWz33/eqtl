# eqtl

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

Map cis/trans expression QTLs from VCF genotypes and a sample×gene phenotype matrix.

## Getting Started

```bash
git clone https://github.com/WWz33/eqtl.git
cd eqtl && make -j

./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff \
  --model lm --mode cis --perm 0 --miss-hand impute -o data/out
```

## Usage

```text
eqtl [options]
```

| Flag | Default | Effect |
|------|---------|--------|
| `-v, --vcf` | required | VCF/BCF genotypes (GT) |
| `-e, --pheno` | required* | phenotype matrix (col1=sample) |
| `-g, --gff` | — | GFF3; omit → genome-wide pairs |
| `-c, --covar` | — | covariates |
| `-k, --grm` | — | relatedness prefix `.grm.id`/`.grm.bin` |
| `--make-grm` | off | write relatedness matrix and exit |
| `-m, --mode` | all | `cis` \| `trans` \| `all` |
| `--model` | lmm | `lm` \| `glm` \| `lmm` \| `glmm` (comma list) |
| `-w, --window` | 1000000 | cis window around TSS (bp) |
| `--pval-cis` | 1e-5 | cis output p threshold |
| `--pval-trans` | 1e-5 | trans/gw output p threshold |
| `--miss-hand` | filter | `filter` \| `impute` missing GT |
| `--fast` | off | share variance/dispersion params per gene |
| `--perm` | 10000 | gene-level permutations (`0`=off) |
| `--seed` | — | permutation seed |
| `--disable-beta-approx` | off | omit beta-approximated p |
| `-o, --out` | eqtl_out | output prefix |
| `-t, --thread` | 1 | threads |

\* not required with `--make-grm`

## Input / Output

| Kind | Notes |
|------|--------|
| Pheno | header = gene ids; col1 = sample; remaining = values |
| Covar | col1 = sample; intercept added if absent |
| GFF | `gene` features; TSS = start (+) or end (−) |
| GRM | `prefix.grm.id` + `prefix.grm.bin` (float32 lower triangle) |
| Pairs | `{out}.{model}.{cis\|trans\|gw}.pairs.tsv` |
| Top | `{out}.{model}.{scope}.top.tsv` |
| Region | `{out}.{model}.{scope}.region.tsv` (`acat_p`, `p_emp`, `p_beta`) |

## Citation

If you use eqtl, cite the methods used in your study design.

## License

MIT
