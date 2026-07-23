#!/usr/bin/env bash
# Build test panel: default 200 samples × 1000 genes from full VCF+GFF.
# Override: N_SAMPLE=200 N_GENE=1000 REG_END=20000000 ./scripts/make_test.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="${OUT_DIR:-$ROOT/data/test}"
SRC_VCF="${SRC_VCF:-/home/ww/eqtl/373_xas_ld5.vcf.gz}"
SRC_GFF="${SRC_GFF:-/home/ww/eqtl/wm82_T2T.gff_agat}"
N_SAMPLE="${N_SAMPLE:-200}"
N_GENE="${N_GENE:-1000}"
# wider region so 1000 genes fit on first contig (bp)
REG_END="${REG_END:-50000000}"
mkdir -p "$DATA"

command -v bcftools >/dev/null || { echo "[E] bcftools required" >&2; exit 1; }
[[ -f "$SRC_VCF" ]] || { echo "[E] missing VCF: $SRC_VCF" >&2; exit 1; }
[[ -f "$SRC_GFF" ]] || { echo "[E] missing GFF: $SRC_GFF" >&2; exit 1; }

bcftools query -l "$SRC_VCF" | head -n "$N_SAMPLE" > "$DATA/test.samples"
ns=$(wc -l < "$DATA/test.samples")
(( ns >= N_SAMPLE )) || echo "[W] only $ns samples available (wanted $N_SAMPLE)" >&2

CONTIG=$(bcftools view -h "$SRC_VCF" | awk -F'[=,>]' '/^##contig/ {print $3; exit}')
REG="${CONTIG}:1-${REG_END}"
echo "[I] region $REG  samples=$ns  genes<=$N_GENE" >&2

bcftools view -S "$DATA/test.samples" -r "$REG" -Ou "$SRC_VCF" \
  | bcftools view -m2 -M2 -v snps -Oz -o "$DATA/test.vcf.gz"
bcftools index -f -t "$DATA/test.vcf.gz"

awk -v c="$CONTIG" -v n="$N_GENE" -v end="$REG_END" 'BEGIN{FS=OFS="\t"}
  $0!~/^#/ && $1==c && $3=="gene" {
    id=""
    n_attr=split($9,a,";")
    for(i=1;i<=n_attr;i++){
      split(a[i],kv,"=")
      if(kv[1]=="ID"){id=kv[2]; gsub(/^gene:/,"",id)}
    }
    if(id=="") next
    tss=($7=="-")?$5:$4
    if(tss>=1 && tss<=end){
      print id, c, tss, $7
      k++
      if(k>=n) exit
    }
  }' "$SRC_GFF" > "$DATA/test.genes.tsv"

ng=$(wc -l < "$DATA/test.genes.tsv")
(( ng > 0 )) || { echo "[E] no genes in region" >&2; exit 1; }
echo "[I] genes written: $ng" >&2

{
  echo '##gff-version 3'
  while read -r id c tss strand; do
    if [[ "$strand" == "-" ]]; then
      s=$((tss-100)); e=$tss
    else
      s=$tss; e=$((tss+100))
    fi
    [[ $s -lt 1 ]] && s=1
    printf '%s\t.\tgene\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s\n' "$c" "$s" "$e" "$strand" "$id" "$id"
  done < "$DATA/test.genes.tsv"
} > "$DATA/test.gff"

python3 - "$DATA" <<'PY'
import random, sys
from pathlib import Path
data = Path(sys.argv[1])
samples = [l.strip() for l in open(data / "test.samples") if l.strip()]
genes = [l.split()[0] for l in open(data / "test.genes.tsv") if l.strip()]
random.seed(1)
# continuous pheno (lm/lmm)
with open(data / "test.pheno.tsv", "w") as f:
    f.write("sample\t" + "\t".join(genes) + "\n")
    for i, s in enumerate(samples):
        row = [s]
        for j, _ in enumerate(genes):
            v = random.gauss(0, 1) + (0.5 * ((i % 3) - 1) if j < 5 else 0.0)
            row.append(f"{v:.6f}")
        f.write("\t".join(row) + "\n")
# counts (glm/glmm)
random.seed(2)
with open(data / "test.counts.tsv", "w") as f:
    f.write("sample\t" + "\t".join(genes) + "\n")
    for i, s in enumerate(samples):
        row = [s]
        for j, _ in enumerate(genes):
            mu = 20 + (3 if j < 5 else 0) + (i % 3)
            row.append(str(max(0, int(random.gauss(mu, 3)))))
        f.write("\t".join(row) + "\n")
# covar
with open(data / "test.covar.tsv", "w") as f:
    f.write("sample_id\tcov1\n")
    for i, s in enumerate(samples):
        f.write(f"{s}\t{i}\n")
print(f"pheno/counts/covar: {len(samples)} x {len(genes)}")
PY

nsnp=$(bcftools view -H "$DATA/test.vcf.gz" | wc -l)
echo "test ready: $DATA"
echo "  samples=$ns  genes=$ng  snps=$nsnp  region=$REG"
echo "  files: test.vcf.gz test.gff test.pheno.tsv test.counts.tsv test.covar.tsv test.samples test.genes.tsv"
