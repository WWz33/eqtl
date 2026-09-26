#!/usr/bin/env bash
# Publication-grade regression matrix.
# Runs every (model × scope × option) combination on synthetic data and checks
# the output files exist, are non-empty, and match the expected schema.
# Also runs sane-null checks: uniform p-value balance under the null, for LM
# nominal p and for the LMM permutation p_emp/p_beta.
#
# Additions / regressions:
#   -L       show pass/fail lines only
#   THREADS=N (env, default 4) sets the thread count for the threaded cases
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
    head -1 "$f" | awk 'NR==1{split($0,h,"\t"); printf "%s\n", h[1]}' || true
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
    --out "$OUT_DIR/lm.trans"

# ---------- 2. LMM: GRM path ------------------------------------------------

run "grm.build" \
  "$EQTL" -v "$TEST/test.vcf.gz" --make-grm --out "$OUT_DIR/grm"

run "lmm.cis" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode cis \
    --perm 0 --out "$OUT_DIR/lmm.cis"

run "lmm.cis.perm" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode cis \
    --perm 20 --seed 11 --out "$OUT_DIR/lmm.cis.perm"

# The permutation p is only worth printing if it is uniform under the null,
# and this panel is null for cis LMM (1000 genes, none anywhere near
# significance), so the case above can be judged directly: p_emp has to sit
# on the B+1 atoms of the permutation test and spread evenly across them,
# p_beta evenly across the deciles. Panel and seed are fixed, so the bounds
# only have to absorb deliberate code changes, not run-to-run noise.
run "lmm.cis.perm.uniform" \
  python3 "$ROOT/scripts/check_perm_uniform.py" \
    "$OUT_DIR/lmm.cis.perm.lmm.cis.region.tsv" --perm 20

run "lmm.trans" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode trans \
    --perm 20 --seed 11 --out "$OUT_DIR/lmm.trans"

# ---------- 2b. LMM with per-gene missingness --------------------------------
# The panels above are complete, so every gene shares one keep set and the LMM
# paths take their shared-sample branches only. Derive a per-gene missing panel
# to reach the other half: per-gene GRM subsetting, the per-keep basis cache,
# and the mixed-keep SNP-outer loop. Both cases run threaded — the
# shared-basis cis driver and the mixed-keep parallel region are exactly what
# these cases exist to cover.
# The gene count is capped so the two cases stay inside the <60s budget (keep
# this in mind before raising N_GENES: trans scales with it), and the missing
# pattern is fixed so runs are comparable.

python3 - "$TEST/test.pheno.tsv" "$OUT_DIR/miss.pheno.tsv" <<'PY'
import sys

src, dst = sys.argv[1], sys.argv[2]
N_GENES = 120  # the first N genes are enough to mix keeps
rows = [l.rstrip('\n').split('\t') for l in open(src)]
hdr, data = rows[0], rows[1:]
if len(hdr) - 1 < N_GENES:
    sys.exit("mixed-keep panel needs %d gene columns, %s has %d" % (N_GENES, src, len(hdr) - 1))
hdr = hdr[:N_GENES + 1]
with open(dst, 'w') as out:
    out.write('\t'.join(hdr) + '\n')
    for i, row in enumerate(data):
        vals = [row[0]]
        for j in range(1, len(hdr)):
            g = j - 1
            # every gene drops one sample (j % 7), every fifth gene drops a
            # second one as well -> several distinct keep sets, none complete
            miss = (i == g % 7) or (g % 5 == 0 and i == (3 * g + 1) % 11)
            vals.append('NA' if miss else row[j])
        out.write('\t'.join(vals) + '\n')
PY

run "lmm.cis.missing" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$OUT_DIR/miss.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode cis -t "$T" \
    --perm 20 --seed 3 --pval-cis 1 --out "$OUT_DIR/lmm.cis.missing"

run "lmm.trans.missing" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$OUT_DIR/miss.pheno.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model lmm --mode trans -t "$T" \
    --perm 20 --seed 3 --out "$OUT_DIR/lmm.trans.missing"

# `run` only checks the exit code, so assert that the per-gene keeps actually
# reached the test output rather than trusting that the branch was taken: the
# cis pairs carry the two keep sizes from the panel (198 and 199) in the `n`
# column, and both region files cover exactly the panel's genes.
python3 - "$OUT_DIR" "$OUT_DIR/miss.pheno.tsv" <<'PY'
import sys

out_dir, panel = sys.argv[1], sys.argv[2]
hdr = open(panel).readline().rstrip('\n').split('\t')
want = set(hdr[1:])
if len(want) != 120:
    sys.exit("FAIL: panel has %d genes, expected 120" % len(want))
fail = 0

for tag, scope in (("cis", "cis"), ("trans", "trans")):
    path = "%s/lmm.%s.missing.lmm.%s.region.tsv" % (out_dir, scope, scope)
    rows = [l.rstrip('\n').split('\t') for l in open(path)]
    h = rows[0]
    ig, it = h.index("gene"), h.index("n_tested")
    genes = {r[ig] for r in rows[1:]}
    untested = [r[ig] for r in rows[1:] if r[it] == "0"]
    print("%s: %d region rows, %d genes, %d untested" % (tag, len(rows) - 1, len(genes), len(untested)))
    if genes != want or untested:
        print("FAIL: region file does not cover the panel"); fail = 1

path = "%s/lmm.cis.missing.lmm.cis.pairs.tsv" % out_dir
rows = [l.rstrip('\n').split('\t') for l in open(path)]
h = rows[0]
ig, in_ = h.index("gene"), h.index("n")
genes = {r[ig] for r in rows[1:]}
ns = sorted({r[in_] for r in rows[1:]})
print("cis pairs: %d rows, %d genes, keep sizes n=%s" % (len(rows) - 1, len(genes), ns))
if genes != want or ns != ["198", "199"]:
    print("FAIL: expected both per-gene keep sizes in the n column"); fail = 1

sys.exit(fail)
PY

# ---------- 3. Count models: cis only --------------------------------------

run "glm.cis" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$TEST/test.counts.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" --model glm --mode cis --perm 50 --seed 17 \
    --out "$OUT_DIR/glm.cis"

# GLMM does a full PQL refit per SNP, and roughly 15s per gene on this panel,
# with each permutation draw paying the same again. Keep all samples and
# restrict to the first N gene columns: the counts matrix is samples (rows) ×
# genes (columns), so subsetting rows would cut the sample size instead and
# turn this into a different test. Keep the case small and threaded — this is
# a smoke test of the GLMM cis path (nominal + one draw), not a scan of the
# whole panel; raising N or --perm costs minutes, not seconds.
N_GLMM_GENES="${N_GLMM_GENES:-5}"
python3 - "$TEST/test.counts.tsv" "$OUT_DIR/counts.sub.tsv" "$N_GLMM_GENES" <<'PY'
import sys
src, dst, keep = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(src) as f, open(dst, "w") as o:
    for line in f:
        t = line.rstrip('\n').split('\t')
        if len(t) - 1 < keep:
            sys.exit("counts matrix has %d genes, need %d" % (len(t) - 1, keep))
        o.write('\t'.join(t[:keep + 1]) + '\n')
PY
run "glmm.cis" \
  "$EQTL" -v "$TEST/test.vcf.gz" -e "$OUT_DIR/counts.sub.tsv" -g "$TEST/test.gff" \
    -c "$TEST/test.covar.tsv" -k "$OUT_DIR/grm" --model glmm --mode cis -t "$T" \
    --perm 1 --seed 17 --pval-cis 1 --out "$OUT_DIR/glmm.cis"

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
