#!/usr/bin/env python3
"""A/B compare two eqtl output trees in three judging modes.

Run both sides with the same --thread value: values are thread-invariant,
row order is not, so rows are matched by their key columns, not by position.

  compare_outputs.py --mode exact A B
  compare_outputs.py --mode tol --rel 1e-12 A B
  compare_outputs.py --mode report A B

A and B are each an output prefix (out.lmm.cis), a directory holding one or
more runs, or a single .tsv. Files are matched by their path relative to the
spec; rows are matched by whichever of --key (default gene,snp) the header
carries.

Judging rules
  exact   every cell must be identical after keying
  tol     integer-valued cells (n, n_tested, pos, ...) must still be
          identical; float cells pass when |a-b| / max(|a|,|b|,1e-300) <= --rel.
          p_emp and p_beta are reported but never judged: they are quantised
          (p_emp moves in steps of 1/(B+1)) and no tolerance band can hold
          them. Non-numeric mismatches, missing rows and row sets that differ
          always fail.
  report  the tol-style drift table only, exit 0 regardless

Exit code: 0 when the mode's verdict passes, 1 otherwise.
"""

import argparse
import glob
import hashlib
import math
import os
import sys

# Permutation-derived columns (quantised p, the beta fit to the permutation
# min-p distribution): report their drift, never judge it. The nominal
# statistics in the same file are still judged as usual.
PERM_COLS = {"p_emp", "p_beta", "beta_shape1", "beta_shape2"}


def collect(spec):
    """basename -> path, from a directory, an output prefix or a single file."""
    if os.path.isdir(spec):
        out = {}
        for root, _, fs in os.walk(spec):
            for f in sorted(fs):
                if f.endswith(".tsv"):
                    p = os.path.join(root, f)
                    out[os.path.relpath(p, spec)] = p
        return out
    if os.path.isfile(spec):
        return {os.path.basename(spec): spec}
    # prefix mode: key on the suffix after the prefix, so two runs whose
    # prefixes differ (out vs out2) still line their files up
    base = os.path.dirname(spec) or "."
    pre = os.path.basename(spec)
    out = {}
    for p in sorted(glob.glob(os.path.join(base, pre + "*.tsv"))):
        out[os.path.basename(p)[len(pre):]] = p
    return out


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f]
    return hdr, rows


def as_number(s):
    """(finite value, is_int_literal) or (None, False) for text.

    NaN and infinities count as text: they must show up as a structural
    difference, not vanish into float comparisons (abs(nan - x) is nan and
    max() quietly drops it, which would turn a real change into a pass).
    """
    try:
        v = float(s)
    except ValueError:
        return None, False
    if not math.isfinite(v):
        return None, False
    return v, ("." not in s and "e" not in s and "E" not in s)


def compare_exact(rows_a, rows_b, hdr):
    if len(rows_a) != len(rows_b):
        return False, ["    row count differs: %d != %d" % (len(rows_a), len(rows_b))]
    bad = 0
    egs = []
    for ra, rb in zip(sorted(rows_a), sorted(rows_b)):
        if ra == rb:
            continue
        bad += 1
        if len(egs) < 3:
            if len(ra) != len(rb):
                egs.append("field count %d != %d" % (len(ra), len(rb)))
            else:
                for c, (x, y) in enumerate(zip(ra, rb)):
                    if x != y:
                        egs.append("%s: %s -> %s" % (hdr[c], x, y))
                        break
    if bad:
        return False, ["    %d/%d rows differ; e.g. %s" % (bad, len(rows_a), "; ".join(egs))]
    return True, ["    all cells identical (row order ignored)"]


def compare_tol(rows_a, rows_b, hdr, keys, rel):
    """Per-column drift over key-aligned rows. Returns (ok, lines, failures)."""
    kcols = [k for k in keys if k in hdr]
    idx = [hdr.index(k) for k in kcols]
    map_a = {}
    map_b = {}
    for rows, m in ((rows_a, map_a), (rows_b, map_b)):
        for r in rows:
            if len(r) != len(hdr):
                return False, ["    malformed row (%d of %d fields)" % (len(r), len(hdr))], 1
            k = tuple(r[i] for i in idx)
            if k in m:
                return False, ["    key %s is not unique — cannot align rows"
                               % ",".join(kcols or keys)], 1
            m[k] = r
    if set(map_a) != set(map_b):
        return False, ["    row sets differ (A=%d B=%d); only-A %s only-B %s"
                       % (len(map_a), len(map_b),
                          sorted(set(map_a) - set(map_b))[:3],
                          sorted(set(map_b) - set(map_a))[:3])], 1

    n_rows = len(map_a)
    lines = []
    failures = 0
    for ci, name in enumerate(hdr):
        if ci in idx:
            continue
        changed = int_bad = text_bad = 0
        mx_abs = mx_rel = 0.0
        for k, ra in map_a.items():
            x, y = ra[ci], map_b[k][ci]
            if x == y:
                continue
            fx, ix = as_number(x)
            fy, iy = as_number(y)
            if fx is None or fy is None:
                text_bad += 1
                continue
            if ix and iy:
                int_bad += 1
                continue
            changed += 1
            d = abs(fx - fy)
            mx_abs = max(mx_abs, d)
            mx_rel = max(mx_rel, d / max(abs(fx), abs(fy), 1e-300))
        if not (changed or int_bad or text_bad):
            continue
        notes = []
        if text_bad:
            notes.append("%d non-numeric" % text_bad)
        if int_bad:
            notes.append("%d integer" % int_bad)
        if changed:
            notes.append("max_abs=%.3g max_rel=%.3g" % (mx_abs, mx_rel))
        structural = text_bad + int_bad > 0
        over = changed > 0 and mx_rel > rel and name not in PERM_COLS
        if name in PERM_COLS:
            tag = "  [reported, not judged]"
        elif structural:
            tag = "  <-- FAIL (structural)"
        elif over:
            tag = "  <-- FAIL (over tolerance)"
        else:
            tag = ""
        if structural or over:
            failures += 1
        lines.append("    %-14s %d/%d rows changed (%s)%s"
                     % (name, changed + int_bad + text_bad, n_rows, ", ".join(notes), tag))
    if not lines:
        lines = ["    all columns identical"]
    return failures == 0, lines, failures


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", help="output prefix, directory or .tsv of side A")
    ap.add_argument("b", help="same for side B")
    ap.add_argument("--mode", choices=("exact", "tol", "report"), default="exact")
    ap.add_argument("--rel", type=float, default=1e-12,
                    help="relative tolerance for --mode tol [1e-12]")
    ap.add_argument("--key", default="gene,snp",
                    help="comma-separated key columns; those present are used [gene,snp]")
    ap.add_argument("--label", default="", help="prefix for the summary line")
    args = ap.parse_args()

    if os.path.isfile(args.a) and os.path.isfile(args.b):
        # two plain files: compare them as a pair whatever they are called
        fa, fb = {"": args.a}, {"": args.b}
    else:
        fa, fb = collect(args.a), collect(args.b)
    if not fa:
        print("no .tsv files matched %s" % args.a)
        return 1
    keys = [k for k in args.key.split(",") if k]

    ok = True
    for rel_path in sorted(fa):
        pa = fa[rel_path]
        name = ("%s %s" % (args.label, os.path.basename(pa))).strip()
        if rel_path not in fb:
            print("  %s: missing on side B" % name)
            ok = False
            continue
        pb = fb[rel_path]
        hdr_a, rows_a = load(pa)
        hdr_b, rows_b = load(pb)
        same_bytes = md5(pa) == md5(pb)
        print("  %-56s bytes=%s" % (name, "same" if same_bytes else "differ"))
        if hdr_a != hdr_b:
            print("    header differs: %s != %s" % (hdr_a, hdr_b))
            ok = False
            continue
        if args.mode == "exact":
            fok, lines = compare_exact(rows_a, rows_b, hdr_a)
        else:
            fok, lines, _ = compare_tol(rows_a, rows_b, hdr_a, keys, args.rel)
        for l in lines:
            print(l)
        if not fok:
            ok = False

    if args.mode == "report":
        print("  verdict: report only (%s)" % ("clean" if ok else "differences found"))
        return 0
    print("  verdict: %s (mode=%s%s)"
          % ("PASS" if ok else "FAIL", args.mode,
             ", rel=%g" % args.rel if args.mode == "tol" else ""))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
