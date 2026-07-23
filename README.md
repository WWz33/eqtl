# eqtl

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

cis/trans eQTL from VCF/BCF or PLINK bed and a sample×gene phenotype matrix.

## Getting Started

```bash
git clone --recurse-submodules https://github.com/WWz33/eqtl.git
cd eqtl && make -j
# make USE_OPENBLAS=0  # pure Eigen if OpenBLAS unavailable

# index genotypes (bcftools; CSI or TBI) when using --vcf
bcftools index -t data/smoke.vcf.gz

# optional: PLINK bed (same panel as GCTA)
plink2 --vcf data/smoke.vcf.gz --make-bed --out data/smoke --allow-extra-chr

# GRM via GCTA
gcta64 --bfile data/smoke --make-grm --out data/smoke_grm

# LMM eQTL from VCF
./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff \
  -k data/smoke_grm --model lmm --mode cis --perm 0 --miss-hand impute \
  -o data/out

# or from PLINK bfile (GCTA/GEMMA-style), with MAF filter
./eqtl -b data/smoke -e data/smoke.pheno.tsv -g data/smoke.gff \
  -k data/smoke_grm --model lmm --mode cis --perm 0 --miss-hand impute --maf 0.05 \
  -o data/out_bed
```

BCF + CSI is preferred for large VCF panels:

```bash
bcftools view -Ob -o panel.bcf panel.vcf.gz
bcftools index panel.bcf
```

## Usage

```text
eqtl [options]
```

| Flag | Default | Effect |
|------|---------|--------|
| `-v, --vcf` | * | VCF/BCF genotypes (GT); exclusive with `--bfile` |
| `-b, --bfile` | * | PLINK bfile prefix |
| `-e, --pheno` | required* | phenotype matrix (col1=sample) |
| `-g, --gff` | — | GFF3 gene features |
| `--gff-id-key` | — | GFF attribute for gene id |
| `-c, --covar` | — | covariates |
| `-k, --grm` | — | relatedness prefix `.grm.id`/`.grm.bin` |
| `--make-grm` | off | write relatedness matrix and exit |
| `-m, --mode` | all | `cis` \| `trans` \| `all` \| `gw` |
| `--model` | lmm | `lm` \| `glm` \| `lmm` \| `glmm` (comma list) |
| `-w, --window` | 1000000 | cis window around TSS (bp) |
| `--pval-cis` | 1e-5 | cis output p threshold |
| `--pval-trans` | 1e-5 | trans/gw output p threshold |
| `--miss-hand` | filter | `filter` \| `impute` missing GT |
| `--max-miss` | 0 | drop SNP if missing fraction > value |
| `--maf` | 0 | min MAF (`0`=off) |
| `--fast` | off | LMM: sparse GRM approx; glm/glmm: fix null phi/sigma2 |
| `--perm` | 0 | gene-level permutations (`0`=off) |
| `--seed` | — | permutation seed |
| `--disable-beta-approx` | off | omit beta-approximated p |
| `-o, --out` | eqtl_out | output prefix |
| `-t, --thread` | 1 | threads |

\* need exactly one of `--vcf` or `--bfile`; pheno not required with `--make-grm`

## Input files

### Genotypes (`-v/--vcf` or `-b/--bfile`)

Exactly one of `--vcf` (VCF/BCF, GT) or `--bfile` (PLINK `.bed`/`.bim`/`.fam` prefix).

Index VCF/BCF with **bcftools**:

```bash
bcftools index -t panel.vcf.gz
# or
bcftools view -Ob -o panel.bcf panel.vcf.gz && bcftools index panel.bcf
```

| `--miss-hand` | `--max-miss` | Effect |
|---------------|--------------|--------|
| `filter` | `0` | any missing → drop SNP |
| `filter` | `m>0` | drop if missing fraction > m; else mean-impute remaining |
| `impute` | `m` | drop if fraction > m; else mean-impute |

### Phenotype (`-e/--pheno`)

TSV, row = sample, columns = genes.

```text
sample	geneA	geneB
S1	1.2	3.4
S2	0.5	2.1
```

| Rule | |
|------|--|
| Header | required |
| Col1 | sample ID |
| Other headers | gene IDs |
| `lm` / `lmm` | continuous values |
| `glm` / `glmm` | non-negative counts |
| Missing | `NA` / `NaN` / `.` → drop that sample for the gene (GCTA-style) |

### Covariates (`-c/--covar`)

```text
sample	cov1	cov2
S1	0	1.2
S2	1	0.3
```

Col1 = sample ID; remaining columns = covariates.

### Annotation (`-g/--gff`)

GFF3 `gene` lines. Gene id from `ID`, else `Name` / `gene_id` (override with `--gff-id-key`).

| Item | Definition |
|------|------------|
| TSS | `+` → start; `−` → end (GFF 1-based) |
| cis interval | `[TSS−W, TSS+W]` (`-w`, bp) |

### Relatedness (`-k/--grm`)

| File | Content |
|------|---------|
| `{prefix}.grm.id` | one sample per line; `FID IID` → IID, or single ID |
| `{prefix}.grm.bin` | float32 lower triangle incl. diagonal; order = `.id` |

`--make-grm` writes the two files and exits.

## Output files

Prefix `-o PREFIX`. Uncompressed TSV per model and scope:

```text
{PREFIX}.{model}.{scope}.pairs.tsv
{PREFIX}.{model}.{scope}.top.tsv
{PREFIX}.{model}.{scope}.region.tsv
```

`{model}` ∈ `lm,glm,lmm,glmm`. `{scope}` ∈ `cis,trans,gw`.

### `{scope}.pairs.tsv`

Rows with `p ≤ --pval-cis` (cis) or `p ≤ --pval-trans` (trans/gw).

| Column | Meaning |
|--------|---------|
| `gene` | gene ID |
| `snp` | variant ID or `chrom:pos:ref:alt` |
| `chrom` | contig |
| `pos` | 1-based position |
| `ref` / `alt` | alleles; beta on alt dosage |
| `maf` | MAF |
| `beta` / `se` / `stat` / `p` | association |
| `r2` | r² |
| `n` | sample size |
| `tss_dist` | `pos − TSS` (bp); `NA` if no TSS |
| `scope` | `cis` / `trans` / `gw` |
| `phi` | NB dispersion (`glm` only) |
| `glm_converged` / `glmm_converged` | 1/0 |

### `{scope}.top.tsv`

Same columns as pairs. At most one row per gene (lowest p among rows that pass the p threshold).

### `{scope}.region.tsv`

| Column | Meaning |
|--------|---------|
| `gene` | gene ID |
| `chrom` | contig |
| `tss` | TSS |
| `n_tested` | SNPs tested |
| `n_sig` | SNPs in pairs |
| `acat_p` | ACAT of SNP p-values |
| `p_emp` | empirical gene p; `NA` if `--perm 0` |
| `p_beta` | beta-approximated gene p; `NA` if off |
| `beta_shape1` / `beta_shape2` | beta fit |

### `--make-grm`

```text
{PREFIX}.grm.id
{PREFIX}.grm.bin
```

## Citation

See repository release notes / paper when available.

## License

MIT
