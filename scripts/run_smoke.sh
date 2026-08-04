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
echo "smoke ok"
