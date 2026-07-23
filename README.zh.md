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

## Input / Output

| 类型 | 说明 |
|------|------|
| 表型 | 表头=基因；第1列=样本 |
| 协变量 | 第1列=样本；无截距则自动加 |
| GFF | `gene`；TSS = start（+）或 end（−） |
| GRM | `prefix.grm.id` + `prefix.grm.bin`（float32 下三角） |
| pairs/top/region | `{out}.{model}.{scope}.*.tsv` |

## Citation

使用时请按研究设计引用相应方法。

## License

MIT
