#!/usr/bin/env python3
# /// script
# dependencies = ["pandas", "numpy", "matplotlib"]
# ///
"""
nquire_depth_spatial.py

Tests two alternatives to "genome-wide tetraploidy" for genes whose nQuire
call is tetraploid (allele balance near 25:75), using only files already in
the repo:

  1. Collapsed paralogs / copy-number gain.  If 25:75 genes are really two
     (or more) diverged copies collapsed onto one reference locus, their read
     depth should be elevated relative to the genome-wide gene median.  A true
     tetraploid has no depth difference between 2n-, 3n- and 4n-called genes.

  2. Segmental aneuploidy / regional duplication.  If calls reflect real
     regional copy-number changes, same-call genes should cluster along
     scaffolds (long runs).  Genome-wide polyploidy or per-gene noise gives
     calls that are spatially independent.  Tested by comparing the fraction
     of adjacent called-gene pairs with the same call against within-scaffold
     permutations.

Inputs (run from the repo root):
  results/nquire2/<sample>.pergene.nquire.tsv.gz
  results/mosdepth/<sample>.regions.bed.gz     chrom start end gene mean_depth
  bed/<strain>.bed                             chrom start end gene

Call rule is the same as nquire_platform_concordance.py (smallest delta-LL;
0/nan missing).  Depth is normalised to the median depth of all genes with
non-zero depth in that sample (norm_depth = 1 means typical single-copy).

Outputs (default results/nquire2/depth_spatial/):
  depth_by_call.tsv        per sample x call: n, median norm_depth, frac >1.5x, frac <0.75x
  spatial_clustering.tsv   per sample: observed vs permuted same-call adjacency, z, p
  pergene_depth_calls.tsv.gz  merged per-gene table (call, margin, relfit, norm_depth, position)
  depth_by_call.pdf        norm_depth distribution per call, per sample

Usage:
  uv run scripts/nquire_depth_spatial.py
  python scripts/nquire_depth_spatial.py --min-margin 5 --nperm 1000
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

PLOIDIES = ["dip", "tri", "tet"]
LONG_SUFFIX = re.compile(r"_(ont|pb)$")


def call_genes(df):
    d = df.copy()
    deltas = d[[f"d_{p}" for p in PLOIDIES]].apply(pd.to_numeric, errors="coerce")
    deltas = deltas.where(deltas != 0)
    called = deltas.notna().any(axis=1)
    filled = deltas.fillna(np.inf)
    d["call"] = np.where(called, filled.idxmin(axis=1).str.replace("d_", "", regex=False), None)
    srt = np.sort(filled.values, axis=1)
    with np.errstate(invalid="ignore"):
        d["margin"] = np.where(called, srt[:, 1] - srt[:, 0], np.nan)
    free = pd.to_numeric(d["free"], errors="coerce")
    d["relfit"] = np.where(called & (free > 0), deltas.min(axis=1) / free, np.nan)
    return d


def load_sample(sample, indir, depthdir, beddir):
    pg = os.path.join(indir, f"{sample}.pergene.nquire.tsv.gz")
    dp = os.path.join(depthdir, f"{sample}.regions.bed.gz")
    if not os.path.exists(dp):
        return None
    d = call_genes(pd.read_csv(pg, sep="\t"))[["gene", "call", "margin", "relfit"]]
    depth = pd.read_csv(dp, sep="\t", header=None, names=["chrom", "start", "end", "gene", "depth"])
    med = depth.loc[depth["depth"] > 0, "depth"].median()
    depth["norm_depth"] = depth["depth"] / med
    m = depth.merge(d, on="gene", how="left")
    # gene order: prefer the strain BED (same coords as regions.bed, kept for safety)
    bedf = os.path.join(beddir, LONG_SUFFIX.sub("", sample) + ".bed")
    if os.path.exists(bedf) and m["start"].isna().any():
        bed = pd.read_csv(bedf, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "gene"])
        m = m.drop(columns=["chrom", "start", "end"]).merge(bed, on="gene", how="left")
    m.insert(0, "sample", sample)
    return m


def depth_by_call(m, min_margin):
    rows = []
    for label, sub in [("all", m), (f"margin>={min_margin:g}", m[m["margin"] >= min_margin])]:
        for c in PLOIDIES:
            x = sub.loc[sub["call"] == c, "norm_depth"].dropna()
            rows.append({
                "sample": m["sample"].iat[0], "subset": label, "call": c, "n": len(x),
                "median_norm_depth": round(x.median(), 3) if len(x) else np.nan,
                "frac_gt1.5x": round(float(np.mean(x > 1.5)), 3) if len(x) else np.nan,
                "frac_lt0.75x": round(float(np.mean(x < 0.75)), 3) if len(x) else np.nan,
            })
    return rows


def same_call_adjacency(calls_by_scaffold):
    same = tot = 0
    for c in calls_by_scaffold:
        if len(c) < 2:
            continue
        same += int(np.sum(c[1:] == c[:-1]))
        tot += len(c) - 1
    return same, tot


def spatial_test(m, min_margin, nperm, rng):
    sub = m[m["call"].notna() & (m["margin"] >= min_margin)].sort_values(["chrom", "start"])
    groups = [g["call"].to_numpy() for _, g in sub.groupby("chrom", sort=False)]
    obs, tot = same_call_adjacency(groups)
    if tot == 0:
        return None
    perm = np.empty(nperm)
    for i in range(nperm):
        perm[i] = same_call_adjacency([rng.permutation(g) for g in groups])[0]
    sd = perm.std(ddof=1)
    # longest run of consecutive tet calls (in called genes) anywhere in the genome
    longest = 0
    for g in groups:
        run = 0
        for c in g:
            run = run + 1 if c == "tet" else 0
            longest = max(longest, run)
    return {
        "sample": m["sample"].iat[0], "genes_used": len(sub), "adjacent_pairs": tot,
        "obs_frac_same": round(obs / tot, 4), "perm_frac_same": round(perm.mean() / tot, 4),
        "z": round((obs - perm.mean()) / sd, 2) if sd > 0 else np.nan,
        "p_upper": (np.sum(perm >= obs) + 1) / (nperm + 1),
        "longest_tet_run": longest,
    }


def plot(merged, outpath, min_margin):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return
    samples = list(merged["sample"].unique())
    ncol = 4
    nrow = int(np.ceil(len(samples) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.6 * nrow), squeeze=False)
    colors = {"dip": "#4C72B0", "tri": "#DD8452", "tet": "#55A868"}
    bins = np.linspace(0, 3, 61)
    for ax, s in zip(axes.flat, samples):
        sub = merged[(merged["sample"] == s) & (merged["margin"] >= min_margin)]
        for c in PLOIDIES:
            x = sub.loc[sub["call"] == c, "norm_depth"].clip(upper=3)
            if len(x):
                ax.hist(x, bins=bins, histtype="step", density=True, color=colors[c],
                        label={"dip": "2n", "tri": "3n", "tet": "4n"}[c])
        ax.axvline(1, color="grey", lw=0.6, ls=":")
        ax.set_title(s, fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in list(axes.flat)[len(samples):]:
        ax.axis("off")
    axes.flat[0].legend(fontsize=7)
    fig.supxlabel("gene depth / median gene depth (clipped at 3)")
    fig.suptitle(f"Per-gene depth by nQuire call (margin >= {min_margin:g})")
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default="results/nquire2")
    ap.add_argument("--depthdir", default="results/mosdepth")
    ap.add_argument("--beddir", default="bed")
    ap.add_argument("--outdir", default=None, help="default: <indir>/depth_spatial")
    ap.add_argument("--min-margin", type=float, default=5.0)
    ap.add_argument("--nperm", type=int, default=500)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    outdir = a.outdir or os.path.join(a.indir, "depth_spatial")
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    merged, dbc, sp = [], [], []
    for f in sorted(glob.glob(os.path.join(a.indir, "*.pergene.nquire.tsv.gz"))):
        s = os.path.basename(f).replace(".pergene.nquire.tsv.gz", "")
        m = load_sample(s, a.indir, a.depthdir, a.beddir)
        if m is None or m["call"].notna().sum() == 0:
            print(f"skipping {s}: no depth file or no called genes")
            continue
        merged.append(m)
        dbc += depth_by_call(m, a.min_margin)
        r = spatial_test(m, a.min_margin, a.nperm, rng)
        if r:
            sp.append(r)

    merged = pd.concat(merged, ignore_index=True)
    merged.to_csv(os.path.join(outdir, "pergene_depth_calls.tsv.gz"), sep="\t", index=False)
    dbc = pd.DataFrame(dbc)
    dbc.to_csv(os.path.join(outdir, "depth_by_call.tsv"), sep="\t", index=False)
    sp = pd.DataFrame(sp)
    sp.to_csv(os.path.join(outdir, "spatial_clustering.tsv"), sep="\t", index=False)
    plot(merged, os.path.join(outdir, "depth_by_call.pdf"), a.min_margin)

    pd.set_option("display.width", 200)
    print(dbc[dbc["subset"] != "all"].to_string(index=False), "\n")
    print(sp.to_string(index=False))
    print(f"\nWrote outputs to {outdir}/")


if __name__ == "__main__":
    main()
