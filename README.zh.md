# eqtl

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

cis/trans eQTL 映射。模型：LM、LMM、NB-GLM、GLMM。

## Getting Started

依赖：C++17、Eigen 3、htslib、OpenBLAS、OpenMP。

```bash
git clone --recurse-submodules https://github.com/WWz33/eqtl.git
cd eqtl && make -j

# GRM
gcta64 --bfile panel --make-grm --out panel_grm

# fission：拆分计数矩阵，估计 PEER 因子
./eqtl fission -e counts.tsv --peer-factors 15 --seed 42 -o fiss

# cis-LMM，PEER 因子作协变量
./eqtl -b panel -e fiss.Y2.tsv -c fiss.factors.tsv -g genes.gff \
  -k panel_grm --model lmm --mode cis -o result

# trans，多模型
./eqtl -b panel -e pheno.tsv -g genes.gff -k panel_grm \
  --model lm,lmm --mode trans --perm 1000 -o result
```

VCF/BCF 输入（`-v`）：建议转 indexed BCF。

```bash
bcftools view -Ob -o panel.bcf panel.vcf.gz
bcftools index panel.bcf
```

## 选项

```
eqtl [options]
eqtl fission [options]
```

| 选项 | 默认 | |
|------|------|-|
| `-b, --bfile` | — | PLINK bed/bim/fam 前缀 |
| `-v, --vcf` | — | VCF/BCF（与 `-b` 互斥） |
| `-e, --pheno` | 必选 | 表型矩阵 |
| `-g, --gff` | — | GFF3 基因注释 |
| `--gff-id-key` | auto | GFF 基因 ID 属性 |
| `-c, --covar` | — | 协变量 |
| `-k, --grm` | — | GRM 前缀（`.grm.id`/`.grm.bin`） |
| `--make-grm` | 关 | 写 GRM 后退出 |
| `-m, --mode` | all | `cis`/`trans`/`gw`/`all` |
| `--model` | lmm | `lm`,`lmm`,`glm`,`glmm`（逗号分隔） |
| `-w, --window` | 1000000 | cis 窗口 ± TSS（bp） |
| `--pval-cis` | 1e-5 | pairs 阈值（cis） |
| `--pval-trans` | 1e-5 | pairs 阈值（trans/gw） |
| `--maf` | 0 | 最小 effect 等位基因频率 |
| `--miss-hand` | impute | `filter`/`impute` |
| `--max-miss` | 0.8 | SNP 缺失率上限 |
| `--fast` | 关 | 稀疏 GRM（LMM）；固定离散参数（GLM/GLMM） |
| `--perm` | 0 | 基因级置换 |
| `--perm-trans-thr` | 1e-5 | trans/gw：stage-2 入选阈值 |
| `--perm-trans-top` | 1000 | trans/gw：stage-2 top-K SNP 数 |
| `--seed` | — | 随机种子 |
| `--disable-beta-approx` | 关 | 跳过 beta 近似 p |
| `-t, --thread` | 1 | 线程数 |
| `-o, --out` | eqtl_out | 输出前缀 |

Fission：

| 选项 | 默认 | |
|------|------|-|
| `--peer-factors` | 10 | 因子数 |
| `--epsilon` | 0.5 | thinning 比例 (0,1) |
| `--fission-max-iter` | 1000 | 最大 VB 迭代 |
| `--fission-tol` | 1e-3 | 收敛阈值 |
| `-c, --covar` | — | PEER 前从 Y1 回归已知协变量 |

## 输入

### 基因型

`-b`：PLINK bed/bim/fam。剂量 = A1 计数。
`-v`：VCF/BCF，GT 字段。区域查询使用 CSI/TBI 索引。

| `--miss-hand` | `--max-miss` | 行为 |
|---------------|--------------|------|
| filter | （任意） | 任一样本缺失则丢 SNP |
| impute | m | 缺失率 > m 则丢；否则均值填补 |

### 表型（`-e`）

```
sample	geneA	geneB
S1	1.2	3.4
S2	0.5	2.1
```

第 1 列 = 样本 ID（列头名称不限）；其余 = 基因。lm/lmm：连续值。glm/glmm：非负整数计数。`NA`/`NaN`/`.` → 该基因排除该样本。

### 协变量（`-c`）

```
sample	cov1	cov2
S1	0	1.2
```

第 1 列 = 样本 ID；其余 = 协变量列。

### GFF（`-g`）

GFF3 `gene` 行。基因 ID 取 `ID`（回退 `Name`/`gene_id`；覆盖：`--gff-id-key`）。TSS：`+`→start，`−`→end。

### GRM（`-k`）

GCTA 格式。`{prefix}.grm.id`（每行一样本），`{prefix}.grm.bin`（float32 下三角含对角线）。兼容 `gcta64 --make-grm` 输出。

## 输出

```
{prefix}.{model}.{scope}.pairs.tsv
{prefix}.{model}.{scope}.top.tsv
{prefix}.{model}.{scope}.region.tsv
```

model ∈ {lm, lmm, glm, glmm}；scope ∈ {cis, trans, gw}。

### pairs

p ≤ 阈值的 SNP–基因对。

| 列 | |
|----|-|
| gene | 基因 ID |
| snp | 变异 ID 或 `chrom:pos:ref:alt` |
| chrom, pos | contig，1-based |
| ref, alt | 等位基因；beta 对应 alt 剂量 |
| af | effect 等位基因频率 |
| beta, se, stat, p | 关联统计量 |
| r2 | 偏 R² |
| n | 样本量 |
| tss_dist | pos − TSS；无注释为 NA |
| scope | cis/trans/gw |
| phi | NB 离散参数（仅 glm） |
| glm_converged, glmm_converged | 1/0 |

### top

列同 pairs。每基因至多一行（过阈值中最小 p 的 SNP）。

### region

每基因一行。

| 列 | |
|----|-|
| gene, chrom, tss | |
| n_tested | 检验 SNP 数 |
| n_sig | pairs 中 SNP 数 |
| acat_p | SNP p 的 ACAT |
| q_bh | 跨基因 BH |
| p_emp | 经验基因 p（`--perm 0` 为 NA） |
| p_beta | beta 近似基因 p |
| beta_shape1, beta_shape2 | beta 拟合参数 |

## Fission

`eqtl fission` 将计数/表达矩阵拆分为条件独立的两半，在一半上估计 PEER 因子，另一半用于下游 eQTL 检验。避免同一数据同时用于混杂估计和关联检验（double-dipping）。

拆分方式：

- 整数计数 → 二项式 thinning：Y1[i,d] ~ Binom(Y[i,d], ε)，Y2 = Y − Y1。由 Poisson 分裂性质 Y1 ⊥ Y2 | λ。
- 非计数（连续）→ Gaussian fission：Y1 = εY + Z，Y2 = (1−ε)Y − Z，Z ~ N(0, σ²ε(1−ε))（逐基因）。Y1 ⊥ Y2 | μ。

PEER（Stegle et al. 2012）：带 ARD 先验的变分贝叶斯因子分析，在 Y1 上运行。若指定 `-c`，先从 Y1 回归已知协变量（Y2 不动）。

输出：`{prefix}.Y1.tsv`、`{prefix}.Y2.tsv`（用作 `-e`）、`{prefix}.factors.tsv`（用作 `-c`）。所有 TSV 首列列头为 `sample_id`。

```bash
./eqtl fission -e counts.tsv --peer-factors 15 --epsilon 0.5 --seed 42 -o fiss
./eqtl -b panel -e fiss.Y2.tsv -c fiss.factors.tsv -g genes.gff \
  -k panel_grm --model lmm --mode cis -o result
```

## 置换

`--perm B`：基因级 min-p 置换。

- **cis**：精确。缓存 cis 窗口全部 SNP 剂量；残差置换 B 次（LM: Freedman–Lane；LMM: GLS 残差）。p_emp = (1 + #{T_perm ≥ T_obs}) / (B+1)。
- **trans/gw**：两阶段（FastQTL 式）。Stage 1：全量 nominal 扫描，每基因保留 top-K SNP。Stage 2：对 obs min-p < `--perm-trans-thr` 的基因，置换后仅重测 top-K。相对全 SNP 置换偏保守。

## License

MIT
