# eqtl

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

cis/trans eQTL from VCF and a sample×gene phenotype matrix.

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
| `-m, --mode` | all | `cis` \| `trans` \| `all` \| `gw` |
| `--model` | lmm | `lm` \| `glm` \| `lmm` \| `glmm` (comma list) |
| `-w, --window` | 1000000 | cis window around TSS (bp) |
| `--pval-cis` | 1e-5 | cis output p threshold |
| `--pval-trans` | 1e-5 | trans/gw output p threshold |
| `--miss-hand` | filter | `filter` \| `impute` missing GT |
| `--max-miss` | 0 | drop SNP if missing fraction > value |
| `--fast` | off | share variance/dispersion params per gene |
| `--perm` | 10000 | gene-level permutations (`0`=off) |
| `--seed` | — | permutation seed |
| `--disable-beta-approx` | off | omit beta-approximated p |
| `-o, --out` | eqtl_out | output prefix |
| `-t, --thread` | 1 | threads |

\* not required with `--make-grm`

## Input files

### Genotypes (`-v/--vcf`)

VCF or BCF. Reads GT only (0/1/2). Ignores DS.

- Multiallelic sites skipped.
- Sample IDs match phenotype / covar / GRM by string.
- Cis: uses TBI/CSI when present; otherwise sequential scan.
- Missing GT:

| `--miss-hand` | `--max-miss` | Effect |
|---------------|--------------|--------|
| `filter` | `0` | any missing → drop SNP |
| `filter` | `m>0` | drop if missing fraction > m; else mean-impute remaining |
| `impute` | `m` | drop if fraction > m; else mean-impute |

### Phenotype (`-e/--pheno`)

TSV: row = sample, columns = genes.

```text
sample	geneA	geneB
S1	1.2	3.4
S2	0.5	2.1
```

- Header required. Col1 = sample ID; other headers = gene IDs.
- `lm`/`lmm`: continuous. `glm`/`glmm`: non-negative counts.
- Gene with any NA: skipped.
- Analysis samples = VCF ∩ pheno (∩ covar if given).

### Covariates (`-c/--covar`, optional)

```text
sample_id	cov1	cov2
S1	0	1.2
S2	1	0.3
```

- Col1 = sample ID. Other columns = covariates.
- Adds intercept if no constant column.
- Samples must be in the analysis set.

### Annotation (`-g/--gff`, optional)

GFF3 `gene` features. Gene id from attribute `ID`, else `Name` / `gene_id` (`--gff-id-key` overrides).

- TSS: `+` → start; `−` → end (GFF 1-based).
- Cis interval: `[TSS−W, TSS+W]` (`-w` bp).
- Chromosome names matched with optional `chr` prefix stripped for equality.
- No GFF → `gw` all-pairs (log line).
- 0 GFF hits for a gene: skip. Multiple hits: warn and skip.

### Relatedness (`-k/--grm`, optional)

| File | Content |
|------|---------|
| `{prefix}.grm.id` | one sample per line; `FID IID` → IID, or single ID |
| `{prefix}.grm.bin` | float32 lower triangle incl. diagonal; order = `.id` |

- `lmm`/`glmm` without `-k`: GRM from VCF on overlap samples.
- `--make-grm`: write the two files and exit.

## Output files

Prefix `-o PREFIX`. Uncompressed TSV. Per model and scope:

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
| `ref` / `alt` | alleles; beta is on alt dosage |
| `maf` | MAF in analysis samples |
| `beta` / `se` / `stat` / `p` | association |
| `r2` | r² |
| `n` | sample size |
| `tss_dist` | `pos − TSS` (bp); `NA` without GFF |
| `scope` | `cis` / `trans` / `gw` |
| `phi` | NB dispersion (`glm` only) |
| `glm_converged` / `glmm_converged` | 1/0 |

### `{scope}.top.tsv`

Same columns as pairs. At most one row per gene: lowest-p SNP if it passes the same p threshold.

### `{scope}.region.tsv`

| Column | Meaning |
|--------|---------|
| `gene` | gene ID |
| `chrom` | contig |
| `tss` | TSS (0 if unknown) |
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

## Sample data

| Path | Contents |
|------|----------|
| `data/smoke.*` | small VCF, pheno, GFF, covar |
| `data/test/` | optional larger panel |

## License

MIT
