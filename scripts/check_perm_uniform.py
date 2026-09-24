#!/usr/bin/env python3
"""Uniformity check for the permutation p-values of a region file.

Under the global null the gene-level p_emp of the permutation test is uniform
on the B+1 atoms {1/(B+1), ..., 1}, and the beta-approximation p_beta is
uniform on [0, 1]. The committed test panel is null for cis LMM (no gene
comes anywhere near significance), so run_tests.sh judges both directly on
the lmm.cis.perm case.

The permutation RNG is seeded, so a given panel, binary and option set always
produce the same file: these thresholds absorb deliberate code changes, not
run-to-run noise. A real defect — a broken shuffle, the wrong draw count, a
stale observed statistic — pushes chi-square into the hundreds. The bounds
are far above what the current binary scores and far below any gross failure.

Usage: check_perm_uniform.py REGION.tsv --perm B
           [--pemp-col p_emp] [--pbeta-col p_beta]
           [--max-chi2-pemp 60] [--max-chi2-pbeta 40]
An empty column name skips that check. Exit 0 when all requested checks pass.
"""

import argparse
import sys


def read_col(path, col):
    """(values, n_skipped) for a numeric column, or (None, 0) when absent."""
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if col not in hdr:
            return None, 0
        i = hdr.index(col)
        vals, skipped = [], 0
        for line in f:
            t = line.rstrip("\n").split("\t")
            try:
                vals.append(float(t[i]))
            except (ValueError, IndexError):
                skipped += 1
    return vals, skipped


def chi2(counts, expect):
    return sum((c - expect) ** 2 / expect for c in counts)


def check_pemp(path, col, perm, limit, min_genes):
    vals, skipped = read_col(path, col)
    if vals is None:
        print("FAIL: no %s column in %s" % (col, path))
        return False
    if len(vals) < min_genes:
        print("FAIL: %s has %d genes, need >= %d" % (col, len(vals), min_genes))
        return False
    # every value has to sit on the B+1 grid of the permutation test
    counts = [0] * (perm + 1)
    off_grid = 0
    for v in vals:
        k = round(v * (perm + 1))
        if not (1 <= k <= perm + 1) or abs(v - k / (perm + 1)) > 1e-9:
            off_grid += 1
            continue
        counts[k - 1] += 1
    if off_grid:
        print("FAIL: %s: %d of %d values are off the B+1 grid for --perm %d"
              % (col, off_grid, len(vals), perm))
        return False
    # the grid has B+1 atoms, and they are not independent: one df goes to the
    # sample-size constraint, so df = B
    c2 = chi2(counts, len(vals) / (perm + 1))
    print("%s: n=%d perm=%d chi2=%.2f df=%d limit=%.0f skipped=%d"
          % (col, len(vals), perm, c2, perm, limit, skipped))
    if c2 > limit:
        print("FAIL: %s is not uniform (chi2 over %.0f)" % (col, limit))
        return False
    return True


def check_pbeta(path, col, limit, min_genes):
    vals, skipped = read_col(path, col)
    if vals is None:
        print("FAIL: no %s column in %s" % (col, path))
        return False
    if len(vals) < min_genes:
        print("FAIL: %s has %d genes, need >= %d" % (col, len(vals), min_genes))
        return False
    bins = [0] * 10
    for v in vals:
        if v < 0.0 or v > 1.0:
            print("FAIL: %s has a value outside [0,1]: %g" % (col, v))
            return False
        bins[int(min(v, 0.9999999) * 10)] += 1
    c2 = chi2(bins, len(vals) / 10)
    print("%s: n=%d deciles=%s chi2=%.2f df=9 limit=%.0f skipped=%d"
          % (col, len(vals), bins, c2, limit, skipped))
    if c2 > limit:
        print("FAIL: %s is not uniform (chi2 over %.0f)" % (col, limit))
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("region", help="region.tsv with the per-gene p-values")
    ap.add_argument("--perm", type=int, required=True,
                    help="the --perm B the case was run with")
    ap.add_argument("--pemp-col", default="p_emp", help="[p_emp]; '' skips")
    ap.add_argument("--pbeta-col", default="p_beta", help="[p_beta]; '' skips")
    ap.add_argument("--max-chi2-pemp", type=float, default=60.0)
    ap.add_argument("--max-chi2-pbeta", type=float, default=40.0)
    ap.add_argument("--min-genes", type=int, default=50)
    args = ap.parse_args()

    ok = True
    if args.pemp_col:
        ok &= check_pemp(args.region, args.pemp_col, args.perm,
                         args.max_chi2_pemp, args.min_genes)
    if args.pbeta_col:
        ok &= check_pbeta(args.region, args.pbeta_col,
                          args.max_chi2_pbeta, args.min_genes)
    print("pass" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
