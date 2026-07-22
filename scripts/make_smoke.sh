#!/usr/bin/env bash
# Build small smoke inputs under data/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/data"
SRC_VCF="${SRC_VCF:-/home/ww/eqtl/373_xas_ld5.vcf.gz}"
SRC_GFF="${SRC_GFF:-/home/ww/eqtl/wm82_T2T.gff_agat}"
N_SAMPLE="${N_SAMPLE:-100}"
N_GENE="${N_GENE:-50}"
mkdir -p "$DATA"

bcftools query -l "$SRC_VCF" | head -n "$N_SAMPLE" > "$DATA/smoke.samples"

CONTIG=$(bcftools view -h "$SRC_VCF" | awk -F'[=,>]' '/^##contig/ {print $3; exit}')
REG="${CONTIG}:1-2000000"
bcftools view -S "$DATA/smoke.samples" -r "$REG" -Ou "$SRC_VCF" \
  | bcftools view -m2 -M2 -v snps -Oz -o "$DATA/smoke.vcf.gz"
bcftools index -f -t "$DATA/smoke.vcf.gz"

awk -v c="$CONTIG" -v n="$N_GENE" 'BEGIN{FS=OFS="\t"}
  $0!~/^#/ && $1==c && $3=="gene" {
    id=""
    n_attr=split($9,a,";")
    for(i=1;i<=n_attr;i++){
      split(a[i],kv,"=")
      if(kv[1]=="ID"){id=kv[2]; gsub(/^gene:/,"",id)}
    }
    if(id=="") next
    tss=($7=="-")?$5:$4
    if(tss>=1 && tss<=2000000){
      print id, c, tss, $7
      k++
      if(k>=n) exit
    }
  }' "$SRC_GFF" > "$DATA/smoke.genes.tsv"

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
  done < "$DATA/smoke.genes.tsv"
} > "$DATA/smoke.gff"

python3 - "$DATA" <<'PY' > "$DATA/smoke.pheno.tsv"
import random, sys
from pathlib import Path
data = Path(sys.argv[1])
samples = [l.strip() for l in open(data/"smoke.samples")]
genes = [l.split()[0] for l in open(data/"smoke.genes.tsv")]
random.seed(1)
print("sample\t" + "\t".join(genes))
for i,s in enumerate(samples):
    row=[s]
    for j,_ in enumerate(genes):
        v = random.gauss(0,1) + (0.8*((i%3)-1) if j==0 else 0)
        row.append(f"{v:.6f}")
    print("\t".join(row))
PY

python3 - "$DATA" <<'PY' > "$DATA/smoke.counts.tsv"
import random, sys
from pathlib import Path
data = Path(sys.argv[1])
samples = [l.strip() for l in open(data/"smoke.samples")]
genes = [l.split()[0] for l in open(data/"smoke.genes.tsv")]
random.seed(2)
print("sample\t" + "\t".join(genes))
for i,s in enumerate(samples):
    row=[s]
    for j,_ in enumerate(genes):
        mu = 20 + (5 if j==0 else 0) + (i%3)
        row.append(str(max(0,int(random.gauss(mu,3)))))
    print("\t".join(row))
PY

{
  echo -e "sample_id\tcov1"
  i=0
  while read -r s; do
    echo -e "${s}\t${i}"
    i=$((i+1))
  done < "$DATA/smoke.samples"
} > "$DATA/smoke.covar.tsv"

echo "smoke ready: $(bcftools view -H "$DATA/smoke.vcf.gz" | wc -l) SNPs, $(wc -l < "$DATA/smoke.samples") samples, $(wc -l < "$DATA/smoke.genes.tsv") genes"
