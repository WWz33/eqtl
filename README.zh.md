# eqtl

<!-- README-I18N:START -->

[English](./README.md) | **中文**

<!-- README-I18N:END -->

从 VCF 基因型与 sample×gene 表型矩阵映射 cis/trans eQTL。

## Getting Started

```bash
git clone https://github.com/WWz33/eqtl.git
cd eqtl && make -j

./scripts/make_smoke.sh
./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff \
  -c data/smoke.covar.tsv --model lm --mode cis --perm 0 --miss-hand impute -o data/out
```

## Usage

```text
eqtl [options]
```

| 选项 | 默认 | 说明 |
|------|------|------|
| `-v, --vcf` | 必选 | VCF/BCF（GT→0/1/2） |
| `-e, --pheno` | 必选* | 表型矩阵（第1列 sample） |
| `-g, --gff` | — | GFF3 取 TSS；省略则全基因组 pair |
| `-c, --covar` | — | 协变量（可选） |
| `-k, --grm` | — | GCTA GRM 前缀 |
| `--make-grm` | 关 | 从 VCF 写 GRM 后退出 |
| `-m, --mode` | all | cis\|trans\|all |
| `--model` | lmm | lm\|glm\|lmm\|glmm（逗号多选） |
| `-w, --window` | 1000000 | TSS 两侧 cis 窗（bp） |
| `--pval-cis` | 1e-5 | cis 写出 p 阈值 |
| `--pval-trans` | 1e-5 | trans/gw 写出 p 阈值 |
| `--miss-hand` | filter | filter\|impute |
| `--fast` | 关 | 每基因共享 nuisance |
| `--perm` | 10000 | 基因级置换（0=关） |
| `--seed` | — | 置换种子 |
| `--disable-beta-approx` | 关 | 仅 p_emp |
| `-o, --out` | eqtl_out | 输出前缀 |
| `-t, --thread` | 1 | 线程数 |

\* 与 `--make-grm` 联用时可不给 pheno

## Input / Output

| 类型 | 说明 |
|------|------|
| 表型 | 表头为基因名；第1列样本 |
| 协变量 | 第1列样本；无截距则自动加 |
| GFF | gene；TSS = +start / −end |
| GRM | `prefix.grm.id` + `prefix.grm.bin` |
| pairs/top/region | `{out}.{model}.{scope}.*.tsv` |

## Citation

开发中。按所用分析方法自行引用。

## License

MIT
