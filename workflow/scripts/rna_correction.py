"""The RNA-seq library correction (config rna_correction:, workflow/rules/rna_correction.smk), on one experiment's
RNA-seq libraries.

Per condition, each library's log2 CPM deviation from the mean of its replicates is fitted as an additive smooth of
log2 transcript length and GC (median regression on natural splines, on the transcripts with a mean count >= min_mean
over the condition; predictors clamped to the fitted range, since natural splines extrapolate linearly). The fit is the
factor a correction would divide out, rescaled so the library keeps its total. The library whose factor varies most
(SD of its log2 over the fitted transcripts) is corrected if that SD is above max_sd; the others are then refitted
against the corrected counts, until none is above it. A corrected library's NumReads are divided by its factor and its
TPM recomputed from them (NumReads / EffectiveLength, scaled to its old TPM total over the same rows); every other file
is copied unchanged.

Usage: rna_correction.py <features.tsv> <max_sd> <min_mean> <df> <scan.tsv> <factors.tsv> <n> <in 1..n> <out 1..n>
                         <condition 1..n>
  features.tsv  Name, length, gc per transcript
  out i         <library>_quant.sf: the library names come from these
  condition i   library i's condition; its replicates are the other libraries of that condition (at least 3 each)
  scan.tsv      one row per fit: condition, round, library, transcripts fitted, SD of the log2 factor, Spearman of the
                deviation with log2 length and with GC before and after the fit, and whether this round corrected that
                library
  factors.tsv   Name + one column per corrected library: the factor its NumReads were divided by
"""
import os
import shutil
import sys
from collections import Counter

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog
from scipy.stats import spearmanr


def ns(x, ref, df):
    # R's splines::ns(ref, df) evaluated at x, in another basis of the same space (ESL eq. 5.4-5.5)
    knots = np.concatenate(([ref.min()], np.quantile(ref, np.arange(1, df) / df), [ref.max()]))

    def d(k):
        return (np.maximum(x - knots[k], 0) ** 3 - np.maximum(x - knots[-1], 0) ** 3) / (knots[-1] - knots[k])

    return np.column_stack([x] + [d(k) - d(len(knots) - 2) for k in range(len(knots) - 2)])


def design(L, G, fit_L, fit_G, df):
    return np.column_stack([np.ones(len(L)), ns(L, fit_L, df), ns(G, fit_G, df)])


def median_regression(X, y):
    # the LP of quantreg::rq(tau = 0.5): min sum(u + v) / 2 s.t. X b + u - v = y, u, v >= 0
    n, p = X.shape
    a_eq = sparse.hstack([sparse.csr_matrix(X), sparse.eye(n), -sparse.eye(n)], format="csr")
    c = np.concatenate([np.zeros(p), np.full(2 * n, 0.5)])
    res = linprog(c, A_eq=a_eq, b_eq=y, bounds=[(None, None)] * p + [(0, None)] * (2 * n), method="highs")
    if res.status != 0:
        sys.exit(f"median regression failed: {res.message}")
    return res.x[:p]


def fit_factor(target, refs, cur, L, G, min_mean, df):
    m = cur[[target] + refs]
    lcpm = np.log2((m + 0.5) / m.sum() * 1e6)
    rows = (m.mean(axis=1) >= min_mean).to_numpy()
    dev = (lcpm[target] - lcpm[refs].mean(axis=1)).to_numpy()[rows]
    fit_L, fit_G = L[rows], G[rows]
    X = design(fit_L, fit_G, fit_L, fit_G, df)
    beta = median_regression(X, dev)
    resid = dev - X @ beta
    clamped = design(np.clip(L, fit_L.min(), fit_L.max()), np.clip(G, fit_G.min(), fit_G.max()), fit_L, fit_G, df)
    f = 2 ** (clamped @ beta)
    counts = cur[target].to_numpy()
    f = f * np.sum(counts / f) / np.sum(counts)  # the library keeps its total
    return f, {"library": target, "transcripts": int(rows.sum()), "sd_log2_factor": np.std(np.log2(f[rows]), ddof=1),
               "spearman_length_before": spearmanr(dev, fit_L)[0], "spearman_length_after": spearmanr(resid, fit_L)[0],
               "spearman_gc_before": spearmanr(dev, fit_G)[0], "spearman_gc_after": spearmanr(resid, fit_G)[0]}


def main():
    args = sys.argv[1:]
    if len(args) < 7:
        sys.exit(__doc__)
    features, max_sd, min_mean, df = args[0], float(args[1]), float(args[2]), int(args[3])
    scan_out, factors_out, n = args[4], args[5], int(args[6])
    if len(args) != 7 + 3 * n:
        sys.exit(__doc__)
    ins, outs, cond = args[7:7 + n], args[7 + n:7 + 2 * n], args[7 + 2 * n:]
    libs = [os.path.basename(o).removesuffix("_quant.sf") for o in outs]
    assert len(set(libs)) == n, libs
    assert min(Counter(cond).values()) >= 3, cond

    feat = pd.read_csv(features, sep="\t", index_col=0)
    q = {lib: pd.read_csv(path, sep="\t") for lib, path in zip(libs, ins)}
    tx = q[libs[0]]["Name"]
    assert all(x["Name"].equals(tx) for x in q.values()) and tx.isin(feat.index).all()
    cur = pd.DataFrame({lib: q[lib]["NumReads"].to_numpy(dtype=float) for lib in libs}, index=tx)
    L = np.log2(feat.loc[tx, "length"].to_numpy(dtype=float))
    G = feat.loc[tx, "gc"].to_numpy(dtype=float)

    scan, factors = [], {}
    for cn in dict.fromkeys(cond):
        members = [lib for lib, c in zip(libs, cond) if c == cn]
        done, rnd = [], 0
        while len(done) < len(members):
            rnd += 1
            fits = {lib: fit_factor(lib, [r for r in members if r != lib], cur, L, G, min_mean, df)
                    for lib in members if lib not in done}
            worst = max(fits, key=lambda lib: fits[lib][1]["sd_log2_factor"])
            fix = fits[worst][1]["sd_log2_factor"] > max_sd
            scan += [{"condition": cn, "round": rnd, **row, "corrected": fix and lib == worst}
                     for lib, (_, row) in fits.items()]
            if not fix:
                break
            cur[worst] = cur[worst] / fits[worst][0]
            factors[worst] = fits[worst][0]
            done.append(worst)
    scan = pd.DataFrame(scan)
    print(scan.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    for lib, path_in, path_out in zip(libs, ins, outs):
        if lib not in factors:
            shutil.copyfile(path_in, path_out)
            continue
        x = q[lib].copy()
        assert (x["EffectiveLength"] > 0).all(), lib
        x["NumReads"] = cur[lib].to_numpy()
        rate = x["NumReads"] / x["EffectiveLength"]
        x["TPM"] = rate / rate.sum() * q[lib]["TPM"].sum()
        x.to_csv(path_out, sep="\t", index=False)
        qs = np.quantile(np.log2(factors[lib]), [0.01, 0.1, 0.5, 0.9, 0.99])
        print(f"{lib} corrected: log2 factor q01 {qs[0]:.2f}, q10 {qs[1]:.2f}, median {qs[2]:.2f}, "
              f"q90 {qs[3]:.2f}, q99 {qs[4]:.2f}")
    scan.to_csv(scan_out, sep="\t", index=False)
    pd.DataFrame({"Name": tx.to_numpy(), **factors}).to_csv(factors_out, sep="\t", index=False)  # none: Name only


if __name__ == "__main__":
    main()
