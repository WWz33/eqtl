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

VCF/BCF。只读 GT（0/1/2）。忽略 DS。

- 多等位位点跳过。
- 样本 ID 与表型/协变量/GRM 字符串一致。
- cis：有 TBI/CSI 则区查，否则顺序扫。
- 缺失 GT：

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

- 须有表头。第1列=样本 ID；其余表头=基因 ID。
- `lm`/`lmm`：连续值。`glm`/`glmm`：非负计数。
- 基因含 NA → 跳过该基因。
- 分析样本 = VCF ∩ pheno（有 covar 再交）。

### 协变量（`-c/--covar`，可选）

```text
sample_id	cov1	cov2
S1	0	1.2
S2	1	0.3
```

- 第1列=样本 ID；其余=协变量。
- 无常数列时加截距。
- 样本须属于分析集。

### 注释（`-g/--gff`，可选）

GFF3 的 `gene`。基因 ID 取属性 `ID`，否则 `Name`/`gene_id`（`--gff-id-key` 可改）。

- TSS：`+`→start，`−`→end（GFF 1-based）。
- cis 区间：`[TSS−W, TSS+W]`（`-w`，bp）。
- 染色体名匹配时去掉可选 `chr` 前缀再比。
- 无 GFF → `gw` 全 pair（日志说明）。
- 0 命中 skip；多命中 warn+skip。

### 亲缘矩阵（`-k/--grm`，可选）

| 文件 | 内容 |
|------|------|
| `{prefix}.grm.id` | 每行一个样本；`FID IID` 用 IID，或单列 ID |
| `{prefix}.grm.bin` | float32 下三角（含对角）；顺序同 `.id` |

- `lmm`/`glmm` 无 `-k`：在 overlap 样本上从 VCF 计算 GRM。
- `--make-grm`：写上述两文件后退出。

## 输出文件

前缀 `-o PREFIX`。明文 TSV。每个 model×scope：

```text
{PREFIX}.{model}.{scope}.pairs.tsv
{PREFIX}.{model}.{scope}.top.tsv
{PREFIX}.{model}.{scope}.region.tsv
```

`{model}` ∈ `lm,glm,lmm,glmm`。`{scope}` ∈ `cis,trans,gw`。

### `{scope}.pairs.tsv`

仅 `p ≤ --pval-cis`（cis）或 `p ≤ --pval-trans`（trans/gw）。

| 列 | 含义 |
|----|------|
| `gene` | 基因 ID |
| `snp` | 变异 ID 或 `chrom:pos:ref:alt` |
| `chrom` | contig |
| `pos` | 1-based 坐标 |
| `ref` / `alt` | 等位基因；beta 对应 alt 剂量 |
| `maf` | 分析样本 MAF |
| `beta` / `se` / `stat` / `p` | 关联 |
| `r2` | r² |
| `n` | 样本量 |
| `tss_dist` | `pos − TSS`（bp）；无 GFF 为 `NA` |
| `scope` | `cis` / `trans` / `gw` |
| `phi` | NB 离散参数（仅 `glm`） |
| `glm_converged` / `glmm_converged` | 1/0 |

### `{scope}.top.tsv`

列同 pairs。每基因至多一行：最小 p 的 SNP，且该 SNP 过同一 p 阈。

### `{scope}.region.tsv`

| 列 | 含义 |
|----|------|
| `gene` | 基因 ID |
| `chrom` | contig |
| `tss` | TSS（未知为 0） |
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

## 示例数据

| 路径 | 内容 |
|------|------|
| `data/smoke.*` | 小规模 VCF、表型、GFF、协变量 |
| `data/test/` | 可选更大面板 |

## License

MIT
