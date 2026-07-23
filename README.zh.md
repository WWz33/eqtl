# eqtl

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

从 VCF 与 sample×gene 表型矩阵做 cis/trans eQTL。

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

| 选项 | 默认 | 说明 |
|------|------|------|
| `-v, --vcf` | 必选 | VCF/BCF 基因型（GT） |
| `-e, --pheno` | 必选* | 表型矩阵（第1列 sample） |
| `-g, --gff` | — | GFF3 gene 特征 |
| `--gff-id-key` | — | GFF 基因 ID 属性名 |
| `-c, --covar` | — | 协变量 |
| `-k, --grm` | — | 亲缘矩阵前缀 `.grm.id`/`.grm.bin` |
| `--make-grm` | 关 | 写亲缘矩阵后退出 |
| `-m, --mode` | all | `cis` \| `trans` \| `all` \| `gw` |
| `--model` | lmm | `lm` \| `glm` \| `lmm` \| `glmm`（逗号多选） |
| `-w, --window` | 1000000 | TSS 两侧 cis 窗（bp） |
| `--pval-cis` | 1e-5 | cis 写出 p 阈值 |
| `--pval-trans` | 1e-5 | trans/gw 写出 p 阈值 |
| `--miss-hand` | filter | 缺失 GT：`filter` \| `impute` |
| `--max-miss` | 0 | 缺失比例超过该值则丢 SNP |
| `--fast` | 关 | 每基因共享方差/离散参数 |
| `--perm` | 10000 | 基因级置换（`0`=关） |
| `--seed` | — | 置换种子 |
| `--disable-beta-approx` | 关 | 不写 beta 近似 p |
| `-o, --out` | eqtl_out | 输出前缀 |
| `-t, --thread` | 1 | 线程数 |

\* 与 `--make-grm` 联用时可不给 pheno

## 输入文件

### 基因型（`-v/--vcf`）

VCF/BCF。仅用 **GT**（0/1/2）。多等位位点不使用。

| `--miss-hand` | `--max-miss` | 效果 |
|---------------|--------------|------|
| `filter` | `0` | 任一缺失 → 丢 SNP |
| `filter` | `m>0` | 缺失比例 > m 则丢；否则对剩余缺失均值填补 |
| `impute` | `m` | 比例 > m 则丢；否则均值填补 |

### 表型（`-e/--pheno`）

TSV：行=样本，列=基因。

```text
sample	geneA	geneB
S1	1.2	3.4
S2	0.5	2.1
```

| 项 | |
|----|--|
| 表头 | 必需 |
| 第1列 | 样本 ID |
| 其余表头 | 基因 ID |
| `lm` / `lmm` | 连续值 |
| `glm` / `glmm` | 非负计数 |

### 协变量（`-c/--covar`）

```text
sample_id	cov1	cov2
S1	0	1.2
S2	1	0.3
```

第1列=样本 ID；其余列=协变量。

### 注释（`-g/--gff`）

GFF3 `gene` 行。基因 ID 取 `ID`，否则 `Name` / `gene_id`（可用 `--gff-id-key`）。

| 项 | 定义 |
|----|------|
| TSS | `+`→start；`−`→end（GFF 1-based） |
| cis 区间 | `[TSS−W, TSS+W]`（`-w`，bp） |

### 亲缘矩阵（`-k/--grm`）

| 文件 | 内容 |
|------|------|
| `{prefix}.grm.id` | 每行一个样本；`FID IID` 用 IID，或单列 ID |
| `{prefix}.grm.bin` | float32 下三角（含对角）；顺序同 `.id` |

`--make-grm` 写上述两文件后退出。

## 输出文件

前缀 `-o PREFIX`。明文 TSV，每个 model×scope：

```text
{PREFIX}.{model}.{scope}.pairs.tsv
{PREFIX}.{model}.{scope}.top.tsv
{PREFIX}.{model}.{scope}.region.tsv
```

`{model}` ∈ `lm,glm,lmm,glmm`。`{scope}` ∈ `cis,trans,gw`。

### `{scope}.pairs.tsv`

`p ≤ --pval-cis`（cis）或 `p ≤ --pval-trans`（trans/gw）的行。

| 列 | 含义 |
|----|------|
| `gene` | 基因 ID |
| `snp` | 变异 ID 或 `chrom:pos:ref:alt` |
| `chrom` | contig |
| `pos` | 1-based 坐标 |
| `ref` / `alt` | 等位基因；beta 对应 alt 剂量 |
| `maf` | MAF |
| `beta` / `se` / `stat` / `p` | 关联 |
| `r2` | r² |
| `n` | 样本量 |
| `tss_dist` | `pos − TSS`（bp）；无 TSS 为 `NA` |
| `scope` | `cis` / `trans` / `gw` |
| `phi` | NB 离散参数（仅 `glm`） |
| `glm_converged` / `glmm_converged` | 1/0 |

### `{scope}.top.tsv`

列同 pairs。每基因至多一行（过 p 阈中最小 p 的 SNP）。

### `{scope}.region.tsv`

| 列 | 含义 |
|----|------|
| `gene` | 基因 ID |
| `chrom` | contig |
| `tss` | TSS |
| `n_tested` | 检验 SNP 数 |
| `n_sig` | pairs 中 SNP 数 |
| `acat_p` | SNP p 的 ACAT |
| `p_emp` | 经验基因 p；`--perm 0` 为 `NA` |
| `p_beta` | beta 近似基因 p；关闭为 `NA` |
| `beta_shape1` / `beta_shape2` | beta 拟合参数 |

### `--make-grm`

```text
{PREFIX}.grm.id
{PREFIX}.grm.bin
```

## Citation

见仓库 release / 论文（如有）。

## License

MIT
