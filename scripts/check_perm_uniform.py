#!/usr/bin/env python3
"""Uniformity check for the permutation p-values of a region file.

Under the global null the gene-level p_emp of the permutation test is uniform
on the B+1 atoms {1/(B+1), ..., 1}, and the beta-approximation p_beta is
uniform on [0, 1]. The committed test panel is null for cis LMM (no gene
comes anywhere near significance), so run_tests.sh judges both directly on
the lmm.cis.perm case.

Both statistics are chi-square against the flat expectation (B+1 cells for
p_emp, 10 deciles for p_beta) and are judged by their upper-tail probability,
so the verdict does not depend on B. The permutation RNG is seeded, so a
given panel, binary and option set always produce the same file: the alpha
only has to absorb deliberate code changes, not run-to-run noise.

Usage: check_perm_uniform.py REGION.tsv --perm B
           [--pemp-col p_emp] [--pbeta-col p_beta] [--alpha 1e-5]
An empty column name skips that check. Exit 0 when all requested checks pass.
"""

import argparse
import math
import sys


def _gammaincc(a, x):
    """Regularised upper incomplete gamma Q(a, x): series for x < a+1, Lentz
    continued fraction otherwise (the usual Numerical Recipes split)."""
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(500):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return 1.0 - total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-15:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


def chi2_sf(x, df):
    """P(chi2_df > x). chi2_df is Gamma(df/2, 2), so the tail is Q(df/2, x/2)."""
    return _gammaincc(0.5 * df, 0.5 * x)


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


def verdict(col, c2, df, n, alpha, extra=""):
    p = chi2_sf(c2, df)
    print("%s: n=%d %schi2=%.2f df=%d p=%.3g" % (col, n, extra, c2, df, p))
    if p < alpha:
        print("FAIL: %s is not uniform (p=%.3g < %.3g)" % (col, p, alpha))
        return False
    return True


def check_pemp(path, col, perm, alpha, min_genes):
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
    return verdict(col, c2, perm, len(vals), alpha,
                   "perm=%d skipped=%d " % (perm, skipped))


def check_pbeta(path, col, alpha, min_genes):
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
    return verdict(col, c2, 9, len(vals), alpha,
                   "deciles=%s skipped=%d " % (bins, skipped))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("region", help="region.tsv with the per-gene p-values")
    ap.add_argument("--perm", type=int, required=True,
                    help="the --perm B the case was run with")
    ap.add_argument("--pemp-col", default="p_emp", help="[p_emp]; '' skips")
    ap.add_argument("--pbeta-col", default="p_beta", help="[p_beta]; '' skips")
    ap.add_argument("--alpha", type=float, default=1e-5,
                    help="fail when the chi-square tail probability drops below this [1e-5]")
    ap.add_argument("--min-genes", type=int, default=50)
    args = ap.parse_args()

    ok = True
    if args.pemp_col:
        ok &= check_pemp(args.region, args.pemp_col, args.perm,
                         args.alpha, args.min_genes)
    if args.pbeta_col:
        ok &= check_pbeta(args.region, args.pbeta_col,
                          args.alpha, args.min_genes)
    print("pass" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
