# eqtl

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

cis/trans eQTL mapping. Models: LM, LMM, NB-GLM, GLMM.

## Getting Started

Dependencies: C++17, Eigen 3, htslib, OpenBLAS, OpenMP.

```bash
git clone --recurse-submodules https://github.com/WWz33/eqtl.git
cd eqtl && make -j

# genotype: VCF → PLINK bed
plink2 --vcf panel.vcf.gz --make-bed --out panel --allow-extra-chr

# GRM
gcta64 --bfile panel --make-grm --out panel_grm

# fission: split counts, estimate PEER factors
./eqtl fission -e counts.tsv --peer-factors 15 --seed 42 -o fiss

# cis-LMM with PEER factors as covariates
./eqtl -b panel -e fiss.Y2.tsv -c fiss.factors.tsv -g genes.gff \
  -k panel_grm --model lmm --mode cis -o result

# trans, multiple models
./eqtl -b panel -e pheno.tsv -g genes.gff -k panel_grm \
  --model lm,lmm --mode trans --perm 1000 -o result
```

VCF/BCF input (`-v`): indexed BCF recommended.

```bash
bcftools view -Ob -o panel.bcf panel.vcf.gz
bcftools index panel.bcf
```

## Options

```
eqtl [options]
eqtl fission [options]
```

| Flag | Default | |
|------|---------|-|
| `-b, --bfile` | — | PLINK bed/bim/fam prefix |
| `-v, --vcf` | — | VCF/BCF (exclusive with `-b`) |
| `-e, --pheno` | required | phenotype matrix |
| `-g, --gff` | — | GFF3 gene annotation |
| `--gff-id-key` | auto | GFF attribute for gene ID |
| `-c, --covar` | — | covariates |
| `-k, --grm` | — | GRM prefix (`.grm.id`/`.grm.bin`) |
| `--make-grm` | off | write GRM and exit |
| `-m, --mode` | all | `cis`/`trans`/`gw`/`all` |
| `--model` | lmm | `lm`,`lmm`,`glm`,`glmm` (comma-separated) |
| `-w, --window` | 1000000 | cis window ± TSS (bp) |
| `--pval-cis` | 1e-5 | pairs threshold (cis) |
| `--pval-trans` | 1e-5 | pairs threshold (trans/gw) |
| `--maf` | 0 | min effect-allele frequency |
| `--miss-hand` | impute | `filter`/`impute` |
| `--max-miss` | 0.8 | SNP missingness cutoff |
| `--fast` | off | sparse GRM (LMM); fixed dispersion (GLM/GLMM) |
| `--perm` | 0 | gene-level permutations |
| `--perm-trans-thr` | 1e-5 | trans/gw: stage-2 entry (obs min-p) |
| `--perm-trans-top` | 1000 | trans/gw: top-K SNPs for stage-2 |
| `--seed` | — | RNG seed |
| `--disable-beta-approx` | off | skip beta-approx p |
| `-t, --thread` | 1 | threads |
| `-o, --out` | eqtl_out | output prefix |

Fission:

| Flag | Default | |
|------|---------|-|
| `--peer-factors` | 10 | number of factors |
| `--epsilon` | 0.5 | thinning fraction (0,1) |
| `--fission-max-iter` | 1000 | max VB iterations |
| `--fission-tol` | 1e-3 | convergence tolerance |
| `-c, --covar` | — | residualize covariates from Y1 before PEER |

## Input

### Genotype

`-b`: PLINK bed/bim/fam. Dosage = A1 count.
`-v`: VCF/BCF, GT field. CSI/TBI index used for region queries.

| `--miss-hand` | `--max-miss` | behavior |
|---------------|--------------|----------|
| filter | (any) | drop SNP if any sample missing |
| impute | m | drop SNP if missingness > m; mean-impute remaining |

### Phenotype (`-e`)

```
sample	geneA	geneB
S1	1.2	3.4
S2	0.5	2.1
```

Col 1 = sample ID (header name ignored); remaining = genes. lm/lmm: continuous. glm/glmm: non-negative counts. `NA`/`NaN`/`.` → sample dropped for that gene.

### Covariates (`-c`)

```
sample	cov1	cov2
S1	0	1.2
```

Col 1 = sample ID; remaining = covariate columns.

### GFF (`-g`)

GFF3 `gene` lines. Gene ID from `ID` (fallback `Name`/`gene_id`; override: `--gff-id-key`). TSS: `+`→start, `−`→end.

### GRM (`-k`)

GCTA format. `{prefix}.grm.id` (one sample/line), `{prefix}.grm.bin` (float32 lower-triangle incl. diagonal). Compatible with `gcta64 --make-grm` output.

## Output

```
{prefix}.{model}.{scope}.pairs.tsv
{prefix}.{model}.{scope}.top.tsv
{prefix}.{model}.{scope}.region.tsv
```

model ∈ {lm, lmm, glm, glmm}; scope ∈ {cis, trans, gw}.

### pairs

SNP–gene pairs with p ≤ threshold.

| Column | |
|--------|-|
| gene | gene ID |
| snp | variant ID or `chrom:pos:ref:alt` |
| chrom, pos | contig, 1-based |
| ref, alt | alleles; beta per alt dosage |
| af | effect-allele frequency |
| beta, se, stat, p | association |
| r2 | partial R² |
| n | sample size |
| tss_dist | pos − TSS; NA if unannotated |
| scope | cis/trans/gw |
| phi | NB dispersion (glm only) |
| glm_converged, glmm_converged | 1/0 |

### top

Same columns. ≤1 row per gene (best SNP passing threshold).

### region

One row per gene.

| Column | |
|--------|-|
| gene, chrom, tss | |
| n_tested | SNPs tested |
| n_sig | SNPs in pairs |
| acat_p | ACAT over SNP p-values |
| q_bh | BH across genes |
| p_emp | empirical gene p (NA if `--perm 0`) |
| p_beta | beta-approx gene p |
| beta_shape1, beta_shape2 | beta fit |

## Fission

`eqtl fission` splits a count/expression matrix into two conditionally independent halves, estimates PEER factors on one half, and outputs the other half for downstream eQTL testing. This avoids using the same data for both confounder estimation and association (double-dipping).

Splitting:

- Integer counts → binomial thinning: Y1[i,d] ~ Binom(Y[i,d], ε), Y2 = Y − Y1. Y1 ⊥ Y2 | λ by Poisson splitting.
- Non-count (continuous) → Gaussian fission: Y1 = εY + Z, Y2 = (1−ε)Y − Z, Z ~ N(0, σ²ε(1−ε)) per gene. Y1 ⊥ Y2 | μ.

PEER (Stegle et al. 2012): variational Bayes factor analysis with ARD priors, run on Y1. If `-c` given, known covariates are residualized from Y1 first (Y2 untouched).

Output: `{prefix}.Y1.tsv`, `{prefix}.Y2.tsv` (use as `-e`), `{prefix}.factors.tsv` (use as `-c`). All TSV files have `sample_id` as the first column header.

```bash
./eqtl fission -e counts.tsv --peer-factors 15 --epsilon 0.5 --seed 42 -o fiss
./eqtl -b panel -e fiss.Y2.tsv -c fiss.factors.tsv -g genes.gff \
  -k panel_grm --model lmm --mode cis -o result
```

## Permutation

`--perm B`: gene-level min-p permutation.

- **cis**: exact. All cis-window SNP dosages cached; residuals shuffled B times (Freedman–Lane for LM; GLS residuals for LMM). p_emp = (1 + #{T_perm ≥ T_obs}) / (B+1).
- **trans/gw**: two-stage (FastQTL-style). Stage 1: full nominal scan, retain top-K SNPs per gene. Stage 2: permute and re-test top-K only for genes with obs min-p < `--perm-trans-thr`. Conservative relative to full-SNP permutation.

## License

MIT
