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

VCF or BCF. Only **GT** is read (hard-call → dosage 0/1/2). **DS is ignored**.

- Prefer indexed files (`.tbi` / `.csi`) for cis region queries.
- Multiallelic sites are skipped.
- Sample IDs must match phenotype (and covar/GRM) by string ID.

Missing GT:

| `--miss-hand` | `--max-miss` | Behavior |
|---------------|--------------|----------|
| `filter` (default) | `0` (default) | any missing → drop SNP |
| `filter` | `m>0` | drop if missing fraction > m; remaining miss mean-imputed |
| `impute` | `m` | drop if fraction > m; else mean-impute |

### Phenotype (`-e/--pheno`)

Tab-separated, **row = sample**, **columns = genes**.

```text
sample	geneA	geneB	...
S1	1.2	3.4	...
S2	0.5	2.1	...
```

- Header row required; col1 = sample ID; remaining headers = gene IDs.
- Values: continuous for `lm`/`lmm`; non-negative counts for `glm`/`glmm`.
- Gene with any NA is skipped. Sample set = VCF ∩ pheno (∩ covar if given).

### Covariates (`-c/--covar`, optional)

```text
sample_id	cov1	cov2
S1	0	1.2
S2	1	0.3
```

- Col1 = sample ID (header name free). Remaining columns = covariates.
- Intercept is added if no constant column is present.
- Samples must overlap the analysis set.

### Annotation (`-g/--gff`, optional)

GFF3. Features with type `gene` (default attribute key: `ID`, else `Name` / `gene_id`; override with `--gff-id-key`).

- **TSS**: strand `+` → start; strand `−` → end (1-based as in GFF).
- Cis window: `[TSS−W, TSS+W]` with `-w` (bp). Chromosome names match VCF (`chr1` vs `1` accepted).
- No GFF → mode forced to genome-wide pairs (`gw`); log notes this.
- Gene ID with 0 GFF hits: skip. Multiple hits: warn + skip.

### Relatedness (`-k/--grm`, optional)

GCTA-style pair:

| File | Content |
|------|---------|
| `{prefix}.grm.id` | one sample per line; `FID IID` → use IID, or single token ID |
| `{prefix}.grm.bin` | float32 lower triangle including diagonal, row-major by sample order in `.id` |

Required for `lmm`/`glmm` if not computed on the fly. Without `-k`, mixed models compute GRM from the VCF on overlap samples. `--make-grm` writes the same format and exits.

## Output files

Prefix `-o PREFIX`. Plain TSV (not compressed). One set per model and scope:

```text
{PREFIX}.{model}.{scope}.pairs.tsv
{PREFIX}.{model}.{scope}.top.tsv
{PREFIX}.{model}.{scope}.region.tsv
```

- `{model}`: `lm` / `glm` / `lmm` / `glmm`
- `{scope}`: `cis` / `trans` / `gw` (depends on `--mode` and GFF)

### `{scope}.pairs.tsv` — SNP–gene tests that pass the p threshold

| Column | Meaning |
|--------|---------|
| `gene` | gene ID |
| `snp` | variant ID (or `chrom:pos:ref:alt`) |
| `chrom` | contig |
| `pos` | 1-based position |
| `ref` / `alt` | alleles; effect on alt dosage |
| `maf` | minor allele frequency in analysis samples |
| `beta` | effect size |
| `se` | standard error |
| `stat` | test statistic |
| `p` | nominal p-value |
| `r2` | variance explained (model-dependent) |
| `n` | sample size |
| `tss_dist` | `pos − TSS` (bp); `NA` if no GFF |
| `scope` | `cis` / `trans` / `gw` |
| `phi` | NB dispersion (`glm` only) |
| `glm_converged` / `glmm_converged` | 1/0 (`glm`/`glmm`) |

Only rows with `p ≤ --pval-cis` (cis) or `p ≤ --pval-trans` (trans/gw) are written.

### `{scope}.top.tsv` — best SNP per gene

Same columns as pairs. One row per gene **if** the best SNP passes the same p threshold; otherwise that gene has no top row.

### `{scope}.region.tsv` — gene-level summary

| Column | Meaning |
|--------|---------|
| `gene` | gene ID |
| `chrom` | contig (from GFF or first SNP) |
| `tss` | TSS (0 if unknown) |
| `n_tested` | SNPs tested |
| `n_sig` | SNPs written to pairs |
| `acat_p` | ACAT combination of SNP p-values |
| `p_emp` | empirical gene p (phenotype permutation); `NA` if `--perm 0` |
| `p_beta` | beta-approximated gene p; `NA` if disabled/off |
| `beta_shape1` / `beta_shape2` | beta fit shapes |

### `--make-grm` outputs

```text
{PREFIX}.grm.id
{PREFIX}.grm.bin
```

Same layout as `-k` input.

## Sample data

| Path | Role |
|------|------|
| `data/smoke.*` | small panel for `make` + first run |
| `data/test/` | larger panel (200 samples × 1000 genes) when present |

## Citation

If you use eqtl, cite the methods used in your study design.

## License

MIT
