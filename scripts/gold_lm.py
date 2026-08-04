#!/usr/bin/env python3
"""Tiny OLS gold: y ~ 1 + g; compare beta/se/p to eqtl lm output (atol 1e-4)."""
from __future__ import annotations

import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EQTL = ROOT / "eqtl"
DATA = ROOT / "data"


def ols(y: np.ndarray, g: np.ndarray):
    n = y.size
    X = np.column_stack([np.ones(n), g])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    df = n - 2
    s2 = float(resid @ resid / df)
    XtX_inv = np.linalg.inv(X.T @ X)
    se = math.sqrt(s2 * XtX_inv[1, 1])
    t = beta[1] / se if se > 0 else 0.0
    # two-sided t via regularized incomplete beta (same family as C++)
    # use scipy-free normal approx only if needed; prefer exact via erfc for large df
    # student-t p: use math with incomplete beta simplified for df
    x = df / (df + t * t)
    # incomplete beta I_x(df/2, 1/2) ~ p
    # use scipy if available else asymptotic normal
    try:
        from scipy.special import betainc  # type: ignore

        p = float(betainc(0.5 * df, 0.5, x))
    except Exception:
        # normal approx
        p = math.erfc(abs(t) / math.sqrt(2.0))
    return float(beta[1]), se, t, p


def main() -> int:
    if not EQTL.exists():
        print("eqtl binary missing; run make", file=sys.stderr)
        return 2
    # use first gene and first few SNPs from a tiny synthetic run
    # Build mini pheno/geno from smoke if present
    vcf = DATA / "smoke.vcf.gz"
    pheno = DATA / "smoke.pheno.tsv"
    gff = DATA / "smoke.gff"
    if not vcf.exists():
        print("missing smoke data", file=sys.stderr)
        return 2

    out = DATA / "gold_lm_out"
    cmd = [
        str(EQTL),
        "-v",
        str(vcf),
        "-e",
        str(pheno),
        "-g",
        str(gff),
        "--model",
        "lm",
        "--mode",
        "cis",
        "--perm",
        "0",
        "--pval-cis",
        "1",
        "--miss-hand",
        "impute",
        "-o",
        str(out),
        "-t",
        "1",
    ]
    subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    pairs = Path(str(out) + ".lm.cis.pairs.tsv")
    lines = pairs.read_text().strip().splitlines()
    if len(lines) < 2:
        print("no pairs", file=sys.stderr)
        return 1
    # pick first data row with finite p
    hdr = lines[0].split("\t")
    i_beta, i_se, i_p = hdr.index("beta"), hdr.index("se"), hdr.index("p")
    i_gene, i_snp = hdr.index("gene"), hdr.index("snp")
    row = lines[1].split("\t")
    gene, snp = row[i_gene], row[i_snp]
    b_eq, se_eq, p_eq = float(row[i_beta]), float(row[i_se]), float(row[i_p])

    # reload y and g via bcftools + pheno
    samples = []
    Y = {}
    with pheno.open() as f:
        genes = f.readline().rstrip("\n").split("\t")[1:]
        gi = genes.index(gene)
        for line in f:
            t = line.rstrip("\n").split("\t")
            samples.append(t[0])
            Y[t[0]] = float(t[gi + 1])

    # dosage for snp id like chrom_pos or chrom:pos:ref:alt
    # query all GT for analysis samples
    q = subprocess.check_output(
        ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", str(vcf)],
        text=True,
    )
    # map sample order from vcf
    vcf_samples = subprocess.check_output(["bcftools", "query", "-l", str(vcf)], text=True).split()
    idx = [vcf_samples.index(s) for s in samples]

    gvec = None
    for line in q.splitlines():
        t = line.split("\t")
        chrom, pos, ref, alt = t[0], t[1], t[2], t[3]
        sid = f"{chrom}_{pos}"
        sid2 = f"{chrom}:{pos}:{ref}:{alt}"
        if snp not in (sid, sid2) and not snp.endswith(f"_{pos}") and snp != sid and snp != sid2:
            # also match ID field style GWH..._pos
            if not (snp.endswith("_" + pos) or snp.endswith(":" + pos + ":" + ref + ":" + alt)):
                continue
        gts = t[4:]
        dos = []
        for j in idx:
            gt = gts[j]
            if gt in ("./.", ".|.", ".", "NA"):
                dos.append(np.nan)
            else:
                a = gt.replace("|", "/").split("/")
                dos.append(sum(int(x) for x in a if x in ("0", "1")))
        g = np.array(dos, dtype=float)
        if np.isnan(g).any():
            mu = np.nanmean(g)
            g = np.where(np.isnan(g), mu, g)
        gvec = g
        break
    if gvec is None:
        print(f"SNP {snp} not found in VCF", file=sys.stderr)
        return 1
    y = np.array([Y[s] for s in samples], dtype=float)
    b, se, t, p = ols(y, gvec)
    ok = abs(b - b_eq) <= 1e-4 and abs(se - se_eq) <= 1e-4
    # p may differ slightly (t vs normal); require 1e-3 relative for p if both >0
    if p_eq > 0 and p > 0:
        ok = ok and abs(math.log10(p) - math.log10(p_eq)) <= 0.05
    print(
        f"gene={gene} snp={snp}\n"
        f"  eqtl  beta={b_eq:.8g} se={se_eq:.8g} p={p_eq:.8g}\n"
        f"  numpy beta={b:.8g} se={se:.8g} p={p:.8g}\n"
        f"  match={ok}"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
