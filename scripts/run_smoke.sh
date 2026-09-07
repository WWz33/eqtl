#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
[[ -f data/smoke.vcf.gz ]] || ./scripts/make_smoke.sh
./eqtl -h >/dev/null
./eqtl -v data/smoke.vcf.gz --make-grm -o data/smoke_grm
./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff -c data/smoke.covar.tsv \
  --model lm --mode cis --perm 0 --pval-cis 1 --miss-hand impute -o data/smoke_lm
test -s data/smoke_lm.lm.cis.pairs.tsv

# lmm + perm on bed-backed input (exercises shared-basis path)
./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff -c data/smoke.covar.tsv \
  -k data/smoke_grm --model lmm --mode cis --perm 20 --seed 7 --pval-cis 1 -o data/smoke_lmm
test -s data/smoke_lmm.lmm.cis.region.tsv

# glm on counts (NB IRLS + score test)
[[ -f data/smoke.counts.tsv ]] && ./eqtl -v data/smoke.vcf.gz -e data/smoke.counts.tsv \
  -g data/smoke.gff -c data/smoke.covar.tsv --model glm --mode cis --pval-cis 1 -o data/smoke_glm

# negative tests: must die
if ./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff \
    --model glm --mode trans -o data/smoke_mustfail 2>/dev/null; then
  echo "[E] expected glm+trans to fail" >&2; exit 1
fi
if ./eqtl -v data/smoke.vcf.gz -e data/smoke.pheno.tsv -g data/smoke.gff \
    --model lm --mode cis --window 1e6 -o data/smoke_mustfail 2>/dev/null; then
  echo "[E] expected --window 1e6 to fail" >&2; exit 1
fi
echo "smoke ok"
