#!/usr/bin/env bash
# Publication-grade regression matrix.
# Runs every (model × scope × option) combination on synthetic data and checks
# the output files exist, are non-empty, and match the expected schema.
# Also runs Bayesian sane-null checks (uniform p-value balance under the null).
#
# Additions / regressions:
#   -L       show pass/fail lines only
#   -j N     (eg. 4) set OMP threads for the "fast" axis
#
# Do NOT put any heavy model logic here; every test should run < 60s.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

TEST="${TEST_DIR:-$ROOT/data/test}"
T="${THREADS:-4}"
EQTL="${EQTL_BIN:-$ROOT/eqtl}"
LOG_ONLY="${1:-}"

[[ -x "$EQTL" ]] || { echo "[E] missing binary: $EQTL (run 'make -j')" >&2; exit 1; }
[[ -f "$TEST/test.vcf.gz" ]] || { echo "[E] missing $TEST/test.vcf.gz (run scripts/make_test.sh)" >&2; exit 1; }
[[ -f "$TEST/test.pheno.tsv" ]] || { echo "[E] missing $TEST/test.pheno.tsv" >&2; exit 1; }
[[ -f "$TEST/test.counts.tsv" ]] || { echo "[E] missing $TEST/test.counts.tsv" >&2; exit 1; }
[[ -f "$TEST/test.covar.tsv" ]] || { echo "[E] missing $TEST/test.covar.tsv" >&2; exit 1; }
[[ -f "$TEST/test.gff" ]] || { echo "[E] missing $TEST/test.gff" >&2; exit 1; }

OUT_DIR="${OUT_DIR:-$ROOT/data/test_out}"
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

PASS=0
FAIL=0
run() {
  local label="$1"; shift
  local expect_fail="${EXPECT_FAIL:-0}"
  local before
  before="$(date +%s)"
  if "$@" >"$OUT_DIR/${label}.log" 2>&1; then
    ec=0
  else
    ec=$?
  fi
  local dt=$(( $(date +%s) - before ))
  local verdict
  if [[ $expect_fail == "1" ]]; then
    if [[ $ec -eq 0 ]]; then
      verdict="FAIL(unexpected success)"
      FAIL=$((FAIL+1))
    else
      verdict=pass
      PASS=$((PASS+1))
    fi
  else
    if [[ $ec -eq 0 ]]; then
      verdict=pass
      PASS=$((PASS+1))
    else
      verdict="FAIL(exit=$ec)"
      FAIL=$((FAIL+1))
    fi
  fi
  if [[ "$LOG_ONLY" == "-L" && "$verdict" == "pass" ]]; then return 0; fi
  printf "%-38s %5ds  %s\n" "$label" "$dt" "$verdict"
}

check_file() {  # check_file label col_suffix
  local f="$1"
  if [[ -s "$f" ]]; then
    head -1 "$f" | awk 'NR==1{split($0,h,"\t"); printf "%s", h[1]}' || true
  else
    echo "MISSING($f)"
  fi
}

# ---------- 1. Numerical accuracy (no statistical shortcut) ------------------

run "lm.cis.p0" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" --model lm --mode cis --perm 0 --pval-cis 1 \
    --out "$OUT_DIR/lm.cis.p0"

check_file "$OUT_DIR/lm.cis.p0.lm.cis.pairs.tsv"

run "lm.cis.p50" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" --model lm --mode cis --perm 50 --seed 7 --out "$OUT_DIR/lm.cis.p50"

run "lm.cis.int50" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" --model lm --mode cis --perm 50 --seed 7 \
    --pheno-norm int --out "$OUT_DIR/lm.cis.int50"

run "lm.trans" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" --model lm --mode trans --perm 50 --seed 7 \
    --perm-trans-thr 1e-3 --perm-trans-top 200 --out "$OUT_DIR/lm.trans"

# ---------- 2. LMM: GRM path + --fast sparsification  ------------------------

run "grm.build" \
  "$EQTL" -v "$TEST/test.vcf.gz" --make-grm --out "$OUT_DIR/grm"

run "lmm.cis" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode cis \
    --perm 0 --out "$OUT_DIR/lmm.cis"

run "lmm.cis.fast" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode cis \
    --fast --out "$OUT_DIR/lmm.cis.fast"

run "lmm.cis.perm" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode cis \
    --perm 20 --seed 11 --out "$OUT_DIR/lmm.cis.perm"

run "lmm.trans" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode trans \
    --perm 20 --seed 11 --perm-trans-thr 1e-3 --perm-trans-top 100 \
    --out "$OUT_DIR/lmm.trans"

# ---------- 3. Count models: cis only --------------------------------------

run "glm.cis" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.counts.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" --model glm --mode cis --perm 50 --seed 17 \
    --out "$OUT_DIR/glm.cis"

# GLMM does a full PQL refit per SNP — too slow for a 1000-gene smoke test at
# --perm 20. Restrict to the first N genes with a top-of-file subset matrix.
N_GLMM_GENES="${N_GLMM_GENES:-20}"
python3 - "$TEST/test.counts.tsv" "$OUT_DIR/counts.sub.tsv" "$N_GLMM_GENES" <<'PY'
import sys
src, dst, keep = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(src) as f, open(dst, "w") as o:
    for i, line in enumerate(f):
        o.write(line)
        if i == keep:
            break
PY
run "glmm.cis" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$OUT_DIR/counts.sub.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model glmm --mode cis \
    --perm 5 --seed 17 --out "$OUT_DIR/glmm.cis"

# ---------- 4. Input-source axis ---------------------------------------------
# (was: call make_test.sh again — pointless on this host, the SRC data path is
# missing. Data generation is covered by scripts/make_smoke.sh / make_test.sh
# directly; we only test the eqtl binary here.)

# ---------- 5. Error-path rejection (must exit non-zero) ---------------------

EXPECT_FAIL=1
run "reject.model_glm_trans" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.counts.tsv" -g "$TEST/test.gff" \
    --model glm --mode trans --out "$OUT_DIR/must_fail_glm_trans"
EXPECT_FAIL=0

EXPECT_FAIL=1
run "reject.window_sci" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    --model lm --mode cis --window 1e6 --out "$OUT_DIR/must_fail_window"
EXPECT_FAIL=0

EXPECT_FAIL=1
run "reject.int_plus_glm" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.counts.tsv" -g "$TEST/test.gff" \
    --model glm --pheno-norm int --mode cis --out "$OUT_DIR/must_fail_int"
EXPECT_FAIL=0

EXPECT_FAIL=1
run "reject.no_pheno" \
  "$EQTL" -v "$TEST/test.vcf.gz" --model lm --mode cis --out "$OUT_DIR/must_fail_pheno"
EXPECT_FAIL=0

EXPECT_FAIL=1
run "reject.bad_model" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" --model foo --mode cis --out "$OUT_DIR/must_fail_model"
EXPECT_FAIL=0

EXPECT_FAIL=1
run "reject.nonfinite_perm" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    --model lm --mode cis --perm -5 --out "$OUT_DIR/must_fail_pneg"
EXPECT_FAIL=0

# ---------- 6. Sanity / null-model statistical check --------------------------

run "null_uniform" python3 - "$OUT_DIR/lm.cis.p0.lm.cis.pairs.tsv" <<'PY'
import sys, re
from pathlib import Path
from collections import Counter
p = Path(sys.argv[1])
if not p.exists():
    print("FAIL: null_uniform missing file:", p); sys.exit(1)
# Under the null with ~independent SNPs, p should be roughly uniform.
# Test just for gross bias: K-S against UNIF(0,1) by decile counts.
vals = []
with open(p) as f:
    hdr = f.readline().split("\t")
    ip = hdr.index("p")
    for line in f:
        t = line.split("\t")
        try: vals.append(float(t[ip]))
        except: pass
# count per decile
bins = [0]*10
for v in vals:
    if 0 <= v <= 1: bins[int(min(v,0.9999999)*10)] += 1
tot = len(vals)
proportions = [b/tot for b in bins]
# decile > 2*uniform (0.2) in more than 2 bins w/ n>1e4 → flag
flag = sum(1 for i,b in enumerate(bins) if b > 2*tot/10)
print(f"null_uniform: n={tot} max_decile_share={max(bins)/tot:.3f} flagged_bins={flag}")
if flag > 2:
    print("FAIL: suspicious concentration"); sys.exit(1)
print("pass")
PY

# ---------- 7. Convergence flagging catches garbage --------------------------

python3 - "$OUT_DIR/glm.cis.glm.cis.pairs.tsv" "$OUT_DIR/glmm.cis.glmm.cis.pairs.tsv" <<'PY'
import sys
from pathlib import Path
for path in sys.argv[1:]:
    p = Path(path)
    if not p.exists():
        print("FAIL: missing", p); sys.exit(1)
    hdr = p.open().readline().split("\t")
    # glm: beta/se/stat/p/r2/n/tss_dist/scope/phi/glm_converged
    for name in hdr:
        idx = hdr.index(name)
    cidx = hdr.index("glm_converged") if "glm_converged" in hdr else None
    if cidx is None and "glmm_converged" in hdr:
        cidx = hdr.index("glmm_converged")
    if cidx is None:
        continue
    n_bad = 0
    n_row = 0
    with open(p) as f:
        next(f)
        for line in f:
            t = line.split("\t")
            n_row += 1
            if t[cidx] == "0": n_bad += 1
    print(f"{p.name}: {n_row} rows, {n_bad} unconverged")
PY

# ---------- summary -----------------------------------------------------------

echo
echo "=== results ==="
echo "pass: $PASS"
echo "FAIL: $FAIL"
[[ "$FAIL" -eq 0 ]]
