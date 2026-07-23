# eqtl

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

从 VCF 基因型与 sample×gene 表型矩阵做 cis/trans eQTL。

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
| `-g, --gff` | — | GFF3；省略则全基因组 pair |
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

VCF/BCF。只读 **GT**（硬分型 → 剂量 0/1/2）。**忽略 DS**。

- cis 区查建议有索引（`.tbi` / `.csi`）。
- 多等位位点跳过。
- 样本 ID 与表型（及 covar/GRM）字符串匹配。

缺失 GT：

| `--miss-hand` | `--max-miss` | 行为 |
|---------------|--------------|------|
| `filter`（默认） | `0`（默认） | 任一缺失 → 丢 SNP |
| `filter` | `m>0` | 缺失比例 > m 则丢；剩余缺失均值填补 |
| `impute` | `m` | 比例 > m 则丢；否则均值填补 |

### 表型（`-e/--pheno`）

TSV，**行=样本，列=基因**：

```text
sample	geneA	geneB	...
S1	1.2	3.4	...
S2	0.5	2.1	...
```

- 必须有表头；第1列样本 ID；其余表头=基因 ID。
- `lm`/`lmm`：连续值；`glm`/`glmm`：非负计数。
- 某基因含 NA → 跳过该基因。分析样本 = VCF ∩ pheno（有 covar 再交）。

### 协变量（`-c/--covar`，可选）

```text
sample_id	cov1	cov2
S1	0	1.2
S2	1	0.3
```

- 第1列样本 ID；其余为协变量。
- 无常数列时自动加截距。
- 样本须与分析集重叠。

### 注释（`-g/--gff`，可选）

GFF3。取 `gene` 特征（属性键默认 `ID`，否则 `Name`/`gene_id`；可用 `--gff-id-key`）。

- **TSS**：`+` 链=start，`−` 链=end（GFF 1-based）。
- cis 窗：`[TSS−W, TSS+W]`（`-w`，bp）。染色体名与 VCF 匹配（`chr1`/`1` 可对上）。
- 无 GFF → 强制全基因组 pair（`gw`），并打日志。
- 0 命中 skip；多命中 warn+skip。

### 亲缘矩阵（`-k/--grm`，可选）

GCTA 风格：

| 文件 | 内容 |
|------|------|
| `{prefix}.grm.id` | 每行一个样本；`FID IID` 用 IID，或单列 ID |
| `{prefix}.grm.bin` | float32 下三角（含对角），顺序与 `.id` 一致 |

`lmm`/`glmm` 无 `-k` 时在 overlap 样本上从 VCF 现算。`--make-grm` 写同格式后退出。

## 输出文件

前缀 `-o PREFIX`。明文 TSV（不压缩）。每个 model × scope 一套：

```text
{PREFIX}.{model}.{scope}.pairs.tsv
{PREFIX}.{model}.{scope}.top.tsv
{PREFIX}.{model}.{scope}.region.tsv
```

- `{model}`：`lm` / `glm` / `lmm` / `glmm`
- `{scope}`：`cis` / `trans` / `gw`（由 `--mode` 与是否有 GFF 决定）

### `{scope}.pairs.tsv` — 过 p 阈的 SNP–基因检验

| 列 | 含义 |
|----|------|
| `gene` | 基因 ID |
| `snp` | 变异 ID（或 `chrom:pos:ref:alt`） |
| `chrom` | contig |
| `pos` | 1-based 坐标 |
| `ref` / `alt` | 等位基因；效应对应 alt 剂量 |
| `maf` | 分析样本次等位频率 |
| `beta` | 效应量 |
| `se` | 标准误 |
| `stat` | 统计量 |
| `p` | 名义 p |
| `r2` | 解释方差（随模型） |
| `n` | 样本量 |
| `tss_dist` | `pos − TSS`（bp）；无 GFF 为 `NA` |
| `scope` | `cis` / `trans` / `gw` |
| `phi` | NB 离散参数（仅 `glm`） |
| `glm_converged` / `glmm_converged` | 是否收敛（`glm`/`glmm`） |

仅 `p ≤ --pval-cis`（cis）或 `p ≤ --pval-trans`（trans/gw）的行写入。

### `{scope}.top.tsv` — 每基因最优 SNP

列同 pairs。仅当该基因最优 SNP 过同一 p 阈时写一行；否则该基因无 top 行。

### `{scope}.region.tsv` — 基因级汇总

| 列 | 含义 |
|----|------|
| `gene` | 基因 ID |
| `chrom` | contig |
| `tss` | TSS（未知为 0） |
| `n_tested` | 检验 SNP 数 |
| `n_sig` | 写入 pairs 的 SNP 数 |
| `acat_p` | SNP p 的 ACAT 合并 |
| `p_emp` | 经验基因 p（表型置换）；`--perm 0` 时为 `NA` |
| `p_beta` | beta 近似基因 p；关闭时为 `NA` |
| `beta_shape1` / `beta_shape2` | beta 拟合参数 |

### `--make-grm` 输出

```text
{PREFIX}.grm.id
{PREFIX}.grm.bin
```

与 `-k` 输入同格式。

## 示例数据

| 路径 | 用途 |
|------|------|
| `data/smoke.*` | 小规模首跑 |
| `data/test/` | 较大面板（200×1000，若已生成） |

## Citation

使用时请按研究设计引用相应方法。

## License

MIT
