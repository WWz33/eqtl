#!/usr/bin/env bash
# Compare eqtl vs quasar on shared models (lm, nb_glm, lmm, p_glmm).
# Inputs: data/test (200 sample x 1000 gene). Uses N_GENE_SUBSET for speed.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TEST="${TEST_DIR:-$ROOT/data/test}"
OUT="${OUT_DIR:-$ROOT/data/compare_quasar}"
EQTL="${EQTL_BIN:-$ROOT/eqtl}"
QUASAR="${QUASAR_BIN:-/home/ww/eqtl/quasar/build/quasar}"
PLINK="${PLINK_BIN:-plink}"
N_GENE_SUBSET="${N_GENE_SUBSET:-50}"
WINDOW="${WINDOW:-1000000}"
MODELS="${MODELS:-lm,nb_glm,lmm}"  # p_glmm optional (slow)

mkdir -p "$OUT"
command -v "$PLINK" >/dev/null || PLINK="$(command -v plink2 || true)"
command -v bcftools >/dev/null
[[ -x "$EQTL" ]] || { echo "[E] missing eqtl binary: $EQTL" >&2; exit 1; }
[[ -x "$QUASAR" ]] || { echo "[E] missing quasar: $QUASAR" >&2; exit 1; }

echo "[I] TEST=$TEST OUT=$OUT N_GENE_SUBSET=$N_GENE_SUBSET MODELS=$MODELS" >&2

# --- subset genes ---
head -n "$N_GENE_SUBSET" "$TEST/test.genes.tsv" > "$OUT/genes.sub.tsv"
mapfile -t GENES < <(cut -f1 "$OUT/genes.sub.tsv")
{
  echo '##gff-version 3'
  while read -r id c tss strand; do
    if [[ "$strand" == "-" ]]; then s=$((tss-100)); e=$tss; else s=$tss; e=$((tss+100)); fi
    [[ $s -lt 1 ]] && s=1
    printf '%s\t.\tgene\t%s\t%s\t.\t%s\t.\tID=%s\n' "$c" "$s" "$e" "$strand" "$id"
  done < "$OUT/genes.sub.tsv"
} > "$OUT/sub.gff"

# pheno subset (eqtl orientation: sample x gene)
python3 - "$TEST/test.pheno.tsv" "$OUT/genes.sub.tsv" "$OUT/eqtl.pheno.tsv" <<'PY'
import sys
from pathlib import Path
pheno, genes_f, out = sys.argv[1:4]
keep = [l.split()[0] for l in open(genes_f) if l.strip()]
with open(pheno) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    genes = hdr[1:]
    idx = [genes.index(g) for g in keep]
    with open(out, "w") as o:
        o.write("sample\t" + "\t".join(keep) + "\n")
        for line in f:
            t = line.rstrip("\n").split("\t")
            o.write(t[0] + "\t" + "\t".join(t[i+1] for i in idx) + "\n")
print("eqtl pheno", len(keep), "genes")
PY

# counts subset for nb_glm / glmm
python3 - "$TEST/test.counts.tsv" "$OUT/genes.sub.tsv" "$OUT/eqtl.counts.tsv" <<'PY'
import sys
pheno, genes_f, out = sys.argv[1:4]
keep = [l.split()[0] for l in open(genes_f) if l.strip()]
with open(pheno) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    genes = hdr[1:]
    idx = [genes.index(g) for g in keep]
    with open(out, "w") as o:
        o.write("sample\t" + "\t".join(keep) + "\n")
        for line in f:
            t = line.rstrip("\n").split("\t")
            o.write(t[0] + "\t" + "\t".join(t[i+1] for i in idx) + "\n")
print("eqtl counts", len(keep), "genes")
PY

cp "$TEST/test.covar.tsv" "$OUT/eqtl.covar.tsv"

# --- PLINK: recode plant contig -> numeric 1 for quasar stoi ---
# keep positions; map contig to 1
bcftools annotate --rename-chrs <(printf '%s\t1\n' "$(bcftools view -h "$TEST/test.vcf.gz" | awk -F'[=,>]' '/^##contig/{print $3; exit}')") \
  -Oz -o "$OUT/test.chr1.vcf.gz" "$TEST/test.vcf.gz"
bcftools index -f -t "$OUT/test.chr1.vcf.gz"

# plink1.9 preferred for classic bed
if command -v plink >/dev/null 2>&1; then
  plink --vcf "$OUT/test.chr1.vcf.gz" --make-bed --out "$OUT/plink" --allow-extra-chr --double-id --keep-allele-order 2>"$OUT/plink.log"
else
  plink2 --vcf "$OUT/test.chr1.vcf.gz" --make-bed --out "$OUT/plink" --allow-extra-chr --double-id 2>"$OUT/plink.log"
fi

# fam IID must match sample ids in pheno/cov/grm
# plink may use FID_IID; force FID=0 IID=sample
python3 - "$OUT/plink.fam" "$TEST/test.samples" <<'PY'
import sys
fam, samples = sys.argv[1:3]
s = [l.strip() for l in open(samples) if l.strip()]
# rewrite fam: FID IID ... use sample list order from VCF fam
rows = open(fam).read().splitlines()
# keep plink order of IIDs but ensure IID is sample name
out = []
for i, line in enumerate(rows):
    t = line.split()
    # plink double-id: FID=IID=sample often
    iid = t[1]
    # if double-id sample_sample, fix
    if "_" in iid and iid.count("_")>=1 and iid not in s:
        # try last token? our samples are plain IDs
        pass
    if iid not in s and t[0] in s:
        iid = t[0]
    out.append(f"0 {iid} 0 0 0 -9")
open(fam, "w").write("\n".join(out) + "\n")
print("fam n=", len(out))
PY

# --- quasar phenotype BED: #chr start end phenotype_id + samples (header sample order = bed columns)
# use gene TSS-1,TSS ; chrom=1
python3 - "$OUT/genes.sub.tsv" "$OUT/eqtl.pheno.tsv" "$OUT/quasar.pheno.bed" <<'PY'
import sys
genes_f, pheno, out = sys.argv[1:4]
genes = []
for line in open(genes_f):
    id, c, tss, strand = line.split()
    genes.append((id, int(tss)))
with open(pheno) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    samples = []
    mat = {g: [] for g in hdr[1:]}
    order = hdr[1:]
    for line in f:
        t = line.rstrip("\n").split("\t")
        samples.append(t[0])
        for i, g in enumerate(order):
            mat[g].append(t[i+1])
with open(out, "w") as o:
    o.write("#chr\tstart\tend\tphenotype_id\t" + "\t".join(samples) + "\n")
    for id, tss in genes:
        start = max(0, tss - 1)
        end = tss
        o.write(f"1\t{start}\t{end}\t{id}\t" + "\t".join(mat[id]) + "\n")
print("quasar bed phenotypes", len(genes), "samples", len(samples))
PY

python3 - "$OUT/genes.sub.tsv" "$OUT/eqtl.counts.tsv" "$OUT/quasar.counts.bed" <<'PY'
import sys
genes_f, pheno, out = sys.argv[1:4]
genes = []
for line in open(genes_f):
    id, c, tss, strand = line.split()
    genes.append((id, int(tss)))
with open(pheno) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    samples = []
    mat = {g: [] for g in hdr[1:]}
    order = hdr[1:]
    for line in f:
        t = line.rstrip("\n").split("\t")
        samples.append(t[0])
        for i, g in enumerate(order):
            mat[g].append(t[i+1])
with open(out, "w") as o:
    o.write("#chr\tstart\tend\tphenotype_id\t" + "\t".join(samples) + "\n")
    for id, tss in genes:
        start = max(0, tss - 1)
        end = tss
        o.write(f"1\t{start}\t{end}\t{id}\t" + "\t".join(mat[id]) + "\n")
print("quasar counts bed", len(genes))
PY

# covar: sample_id + intercept + cov1 (quasar adds intercept in v1.1 - still provide cov1)
# use same order as pheno header samples
python3 - "$OUT/eqtl.pheno.tsv" "$OUT/eqtl.covar.tsv" "$OUT/quasar.cov.tsv" <<'PY'
import sys
pheno, cov, out = sys.argv[1:4]
samples = []
with open(pheno) as f:
    f.readline()
    for line in f:
        samples.append(line.split("\t", 1)[0])
cmap = {}
with open(cov) as f:
    f.readline()
    for line in f:
        t = line.split()
        cmap[t[0]] = t[1]
with open(out, "w") as o:
    o.write("sample_id\tcov1\n")
    for s in samples:
        o.write(f"{s}\t{cmap[s]}\n")
print("cov", len(samples))
PY

# GRM for both: eqtl make-grm then convert to quasar dense TSV (sample order = pheno)
"$EQTL" -v "$TEST/test.vcf.gz" -e "$OUT/eqtl.pheno.tsv" --make-grm --miss-hand impute -o "$OUT/eqtl_grm" \
  >"$OUT/eqtl_grm.log" 2>&1
python3 - "$OUT/eqtl_grm.grm.id" "$OUT/eqtl_grm.grm.bin" "$OUT/eqtl.pheno.tsv" "$OUT/quasar.grm.tsv" <<'PY'
import sys, struct
import numpy as np
id_path, bin_path, pheno, out = sys.argv[1:5]
ids = []
for line in open(id_path):
    t = line.split()
    ids.append(t[-1])
n = len(ids)
ntri = n * (n + 1) // 2
raw = open(bin_path, "rb").read()
vals = struct.unpack(f"<{ntri}f", raw[: ntri * 4])
K = np.zeros((n, n), dtype=float)
k = 0
for i in range(n):
    for j in range(i + 1):
        K[i, j] = K[j, i] = vals[k]
        k += 1
# reorder to pheno sample order
samples = []
with open(pheno) as f:
    f.readline()
    for line in f:
        samples.append(line.split("\t", 1)[0])
ix = [ids.index(s) for s in samples]
K2 = K[np.ix_(ix, ix)]
with open(out, "w") as o:
    o.write("sample_id\t" + "\t".join(samples) + "\n")
    for i, s in enumerate(samples):
        o.write(s + "\t" + "\t".join(f"{K2[i,j]:.10g}" for j in range(len(samples))) + "\n")
print("grm", K2.shape)
PY

# --- run eqtl models ---
run_eqtl() {
  local model="$1" pheno="$2" extra=()
  case "$model" in
    lmm|glmm) extra+=(-k "$OUT/eqtl_grm") ;;
  esac
  echo "[I] eqtl model=$model" >&2
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$pheno" -g "$OUT/sub.gff" -c "$OUT/eqtl.covar.tsv" \
    --model "$model" --mode cis --window "$WINDOW" --perm 0 --pval-cis 1 --miss-hand impute --max-miss 1 \
    -o "$OUT/eqtl_${model}" -t "${THREADS:-4}" "${extra[@]}" \
    >"$OUT/eqtl_${model}.log" 2>&1 || { echo "[E] eqtl $model failed"; tail -20 "$OUT/eqtl_${model}.log"; return 1; }
}

# --- run quasar models ---
run_quasar() {
  local model="$1" bed="$2" extra=()
  case "$model" in
    lmm|p_glmm|nb_glmm) extra+=(-g "$OUT/quasar.grm.tsv") ;;
  esac
  echo "[I] quasar model=$model" >&2
  "$QUASAR" -p "$OUT/plink" -b "$bed" -c "$OUT/quasar.cov.tsv" \
    --model "$model" --mode cis -w "$WINDOW" -o "$OUT/quasar_${model}" --verbose \
    "${extra[@]}" >"$OUT/quasar_${model}.log" 2>&1 \
    || { echo "[E] quasar $model failed"; tail -40 "$OUT/quasar_${model}.log"; return 1; }
}

IFS=',' read -ra ML <<< "$MODELS"
for m in "${ML[@]}"; do
  m=$(echo "$m" | tr -d ' ')
  case "$m" in
    lm)
      run_eqtl lm "$OUT/eqtl.pheno.tsv"
      run_quasar lm "$OUT/quasar.pheno.bed"
      ;;
    nb_glm|glm)
      run_eqtl glm "$OUT/eqtl.counts.tsv"
      run_quasar nb_glm "$OUT/quasar.counts.bed"
      ;;
    lmm)
      run_eqtl lmm "$OUT/eqtl.pheno.tsv"
      run_quasar lmm "$OUT/quasar.pheno.bed"
      ;;
    p_glmm|glmm)
      run_eqtl glmm "$OUT/eqtl.counts.tsv"
      run_quasar p_glmm "$OUT/quasar.counts.bed"
      ;;
    *) echo "[W] skip unknown model $m" >&2 ;;
  esac
done

# --- join & report ---
python3 - "$OUT" <<'PY'
import sys, math, re
from pathlib import Path
from collections import defaultdict
out = Path(sys.argv[1])

def read_eqtl_pairs(path):
    # gene snp chrom pos ref alt maf beta se stat p r2 n model mode ...
    rows = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        i = {h: k for k, h in enumerate(hdr)}
        for line in f:
            t = line.rstrip("\n").split("\t")
            gene, snp = t[i["gene"]], t[i["snp"]]
            # normalize snp key to chrom:pos if possible
            pos = t[i["pos"]]
            chrom = t[i["chrom"]]
            key = (gene, str(pos))  # pos-based join (chrom remapped in quasar)
            rows[key] = {
                "beta": float(t[i["beta"]]),
                "se": float(t[i["se"]]),
                "p": float(t[i["p"]]),
                "snp": snp,
                "chrom": chrom,
            }
    return rows

def read_quasar_variant(path):
    # find first *variant* file
    rows = {}
    p = Path(path)
    if not p.exists():
        # quasar naming varies
        cands = list(p.parent.glob(p.name + "*"))
        files = [x for x in cands if x.is_file() and "variant" in x.name.lower()]
        if not files:
            files = [x for x in cands if x.is_file() and x.suffix in (".tsv", ".txt", "")]
        if not files:
            return rows, None
        p = files[0]
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        # flexible names
        def find(*names):
            for n in names:
                if n in hdr:
                    return hdr.index(n)
            # case insensitive
            low = {h.lower(): i for i, h in enumerate(hdr)}
            for n in names:
                if n.lower() in low:
                    return low[n.lower()]
            return None
        ig = find("phenotype_id", "gene", "pheno", "phenotype")
        ip = find("pos", "position", "snp_pos")
        ib = find("beta", "slope", "effect")
        ise = find("se", "std.error", "stderr")
        ipval = find("p", "pval", "p-value", "p_value", "pvalue")
        isnp = find("variant_id", "snp", "rsid", "id")
        if ig is None or ip is None or ib is None or ipval is None:
            return rows, (p, hdr)
        for line in f:
            t = line.rstrip("\n").split("\t")
            gene = t[ig]
            pos = t[ip]
            key = (gene, str(pos))
            se = float(t[ise]) if ise is not None else float("nan")
            rows[key] = {
                "beta": float(t[ib]),
                "se": se,
                "p": float(t[ipval]),
                "snp": t[isnp] if isnp is not None else pos,
            }
    return rows, p

def corr(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0 or dy == 0:
        return float("nan")
    return num / (dx * dy)

def mad(xs):
    if not xs:
        return float("nan")
    xs = sorted(xs)
    mid = xs[len(xs) // 2]
    d = sorted(abs(x - mid) for x in xs)
    return d[len(d) // 2]

pairs = [
    ("lm", "eqtl_lm.lm.cis.pairs.tsv", "quasar_lm"),
    ("nb_glm", "eqtl_glm.glm.cis.pairs.tsv", "quasar_nb_glm"),
    ("lmm", "eqtl_lmm.lmm.cis.pairs.tsv", "quasar_lmm"),
    ("p_glmm", "eqtl_glmm.glmm.cis.pairs.tsv", "quasar_p_glmm"),
]

report = out / "compare_summary.tsv"
with open(report, "w") as rep:
    rep.write("model\tn_eqtl\tn_quasar\tn_join\tcorr_beta\tcorr_log10p\tmedian_abs_dbeta\tmedian_abs_dlog10p\tfrac_sign_agree\tnote\n")
    for model, efile, qprefix in pairs:
        ep = out / efile
        if not ep.exists():
            continue
        er = read_eqtl_pairs(ep)
        # quasar output files
        qfiles = list(out.glob(qprefix + "*"))
        qvar = None
        for f in qfiles:
            if f.is_file() and "variant" in f.name.lower():
                qvar = f
                break
        if qvar is None:
            # any non-log file
            for f in qfiles:
                if f.is_file() and f.suffix in (".tsv", ".txt") and "log" not in f.name:
                    qvar = f
                    break
        if qvar is None:
            rep.write(f"{model}\t{len(er)}\t0\t0\tNA\tNA\tNA\tNA\tNA\tno_quasar_out\n")
            print(f"[{model}] eqtl={len(er)} quasar=MISSING")
            continue
        qr, meta = read_quasar_variant(qvar)
        if not qr:
            # dump header for debug
            with open(qvar) as f:
                hdr = f.readline().strip()
            rep.write(f"{model}\t{len(er)}\t?\t0\tNA\tNA\tNA\tNA\tNA\tparse_fail:{qvar.name}:{hdr[:80]}\n")
            print(f"[{model}] parse fail {qvar} hdr={hdr[:120]}")
            continue
        keys = sorted(set(er) & set(qr))
        betas_e, betas_q, lp_e, lp_q, dbeta, dlp, signs = [], [], [], [], [], [], 0
        for k in keys:
            be, bq = er[k]["beta"], qr[k]["beta"]
            pe, pq = max(er[k]["p"], 1e-300), max(qr[k]["p"], 1e-300)
            betas_e.append(be)
            betas_q.append(bq)
            le, lq = -math.log10(pe), -math.log10(pq)
            lp_e.append(le)
            lp_q.append(lq)
            dbeta.append(abs(be - bq))
            dlp.append(abs(le - lq))
            if (be >= 0) == (bq >= 0):
                signs += 1
        n = len(keys)
        cb = corr(betas_e, betas_q)
        cl = corr(lp_e, lp_q)
        frac = signs / n if n else float("nan")
        rep.write(
            f"{model}\t{len(er)}\t{len(qr)}\t{n}\t{cb:.6g}\t{cl:.6g}\t{mad(dbeta):.6g}\t{mad(dlp):.6g}\t{frac:.6g}\t{qvar.name}\n"
        )
        print(
            f"[{model}] join={n} corr_beta={cb:.4f} corr_-log10p={cl:.4f} "
            f"med|dβ|={mad(dbeta):.4g} med|dlogp|={mad(dlp):.4g} sign_agree={frac:.3f} file={qvar.name}"
        )

print("wrote", report)
# list quasar outputs for user
print("quasar outputs:")
for f in sorted(out.glob("quasar_*")):
    if f.is_file():
        print(" ", f.name, f.stat().st_size)
PY

echo "[I] done. summary: $OUT/compare_summary.tsv" >&2
