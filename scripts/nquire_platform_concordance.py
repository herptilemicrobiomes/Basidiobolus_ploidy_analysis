#!/usr/bin/env python3
"""
nquire_platform_concordance.py

Summarise nQuire per-gene and whole-genome ploidy calls, and test whether
calls agree between short-read and long-read (ONT / PacBio) alignments of
the same strain.

Inputs (default: results/nquire2/):
  <sample>.pergene.nquire.tsv.gz   gene, free, dip, tri, tet, d_dip, d_tri, d_tet
  <sample>.nquire.tsv              whole-genome nQuire lrdmodel output
Sample names ending in _ont or _pb are treated as long-read runs and paired
with the short-read run of the same strain (name without the suffix).

Call rule (matches pipeline/plot/plot_pergene_ploidy.R): winner = smallest
delta-log-likelihood among d_dip, d_tri, d_tet; 0 / nan / -nan are missing;
a gene with all three missing gets no call.

Extra per-gene diagnostics:
  margin  = second-smallest delta minus smallest (how decisive the call is)
  relfit  = best delta / free log-likelihood (0 = fixed model as good as
            the free mixture; near 1 = no fixed ploidy model fits)

Outputs (default: results/nquire2/concordance/):
  pergene_call_summary.tsv      per-sample call fractions, all and confident
  wholegenome_summary.tsv       per-sample whole-genome call and relative fit
  platform_concordance.tsv      per-strain agreement and Cohen's kappa
  crosstab_<strain>.tsv         short (rows) x long (cols) call tables
  platform_shift.pdf            short vs long call composition per strain

Usage:
  python scripts/nquire_platform_concordance.py
  python scripts/nquire_platform_concordance.py --indir results/nquire2 \
      --outdir results/nquire2/concordance --min-margin 5
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
    """Add call, margin and relfit columns to an nQuire per-gene/whole-genome table."""
    d = df.copy()
    deltas = d[[f"d_{p}" for p in PLOIDIES]].apply(pd.to_numeric, errors="coerce")
    deltas = deltas.where(deltas != 0)          # 0 means unestimable
    called = deltas.notna().any(axis=1)
    filled = deltas.fillna(np.inf)
    d["call"] = np.where(called, filled.idxmin(axis=1).str.replace("d_", "", regex=False), None)
    srt = np.sort(filled.values, axis=1)
    with np.errstate(invalid="ignore"):
        d["margin"] = np.where(called, srt[:, 1] - srt[:, 0], np.nan)
    free = pd.to_numeric(d["free"], errors="coerce")
    d["relfit"] = np.where(called & (free > 0), deltas.min(axis=1) / free, np.nan)
    return d


def cohen_kappa(a, b):
    """Unweighted Cohen's kappa for two label vectors (no sklearn dependency)."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) == 0:
        return np.nan
    labels = np.union1d(a, b)
    po = np.mean(a == b)
    pe = sum(np.mean(a == l) * np.mean(b == l) for l in labels)
    return np.nan if pe == 1 else (po - pe) / (1 - pe)


def composition(calls):
    c = pd.Series(calls).dropna().value_counts()
    n = int(c.sum())
    return n, {p: (c.get(p, 0) / n if n else np.nan) for p in PLOIDIES}


def load_pergene(indir):
    out = {}
    for f in sorted(glob.glob(os.path.join(indir, "*.pergene.nquire.tsv.gz"))):
        name = os.path.basename(f).replace(".pergene.nquire.tsv.gz", "")
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            continue
        out[name] = call_genes(df).set_index("gene")
    return out


def load_wholegenome(indir):
    rows = []
    for f in sorted(glob.glob(os.path.join(indir, "*.nquire.tsv"))):
        base = os.path.basename(f)
        if ".perchrom." in base or ".pergene." in base:
            continue
        name = base.replace(".nquire.tsv", "")
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            rows.append({"sample": name, "status": "empty"})
            continue
        r = call_genes(df).iloc[0]
        rows.append({
            "sample": name,
            "status": "ok" if r["call"] is not None else "no_estimate",
            "platform": platform_of(name),
            "free": r.get("free"),
            **{f"d_{p}": r.get(f"d_{p}") for p in PLOIDIES},
            "call": r["call"],
            "margin": r["margin"],
            "relfit": r["relfit"],
        })
    return pd.DataFrame(rows)


def platform_of(name):
    m = LONG_SUFFIX.search(name)
    return {"ont": "ONT", "pb": "PacBio"}[m.group(1)] if m else "short"


def summarise_pergene(res, min_margin):
    rows = []
    for name, d in res.items():
        n_all, comp_all = composition(d["call"])
        conf = d.loc[d["margin"] >= min_margin, "call"]
        n_conf, comp_conf = composition(conf)
        m = d.loc[d["call"].notna(), "margin"]
        rows.append({
            "sample": name,
            "platform": platform_of(name),
            "total_genes": len(d),
            "called": n_all,
            **{f"frac_{p}": round(comp_all[p], 3) for p in PLOIDIES},
            f"confident_n(margin>={min_margin:g})": n_conf,
            **{f"conf_frac_{p}": round(comp_conf[p], 3) for p in PLOIDIES},
            "median_margin": round(m.median(), 2) if n_all else np.nan,
            "frac_margin_lt2": round(float(np.mean(m < 2)), 3) if n_all else np.nan,
            "median_relfit": round(d.loc[d["call"].notna(), "relfit"].median(), 3) if n_all else np.nan,
        })
    return pd.DataFrame(rows)


def concordance(res, min_margin, outdir):
    rows = []
    for long_name in sorted(n for n in res if LONG_SUFFIX.search(n)):
        short_name = LONG_SUFFIX.sub("", long_name)
        if short_name not in res:
            continue
        s, l = res[short_name], res[long_name]
        j = s[["call", "margin"]].join(l[["call", "margin"]], lsuffix="_short", rsuffix="_long", how="inner")
        j = j.dropna(subset=["call_short", "call_long"])
        if j.empty:
            print(f"skipping {short_name}: no genes called in both {short_name} and {long_name}")
            continue
        jc = j[(j["margin_short"] >= min_margin) & (j["margin_long"] >= min_margin)]
        ct = pd.crosstab(j["call_short"], j["call_long"]).reindex(index=PLOIDIES, columns=PLOIDIES, fill_value=0)
        ct.index.name, ct.columns.name = "short_call", f"{platform_of(long_name)}_call"
        ct.to_csv(os.path.join(outdir, f"crosstab_{short_name}.tsv"), sep="\t")
        _, cs = composition(j["call_short"])
        _, cl = composition(j["call_long"])
        rows.append({
            "strain": short_name,
            "long_platform": platform_of(long_name),
            "genes_both_called": len(j),
            "agreement": round(float(np.mean(j["call_short"] == j["call_long"])), 3) if len(j) else np.nan,
            "kappa": round(cohen_kappa(j["call_short"], j["call_long"]), 3),
            "genes_both_confident": len(jc),
            "agreement_confident": round(float(np.mean(jc["call_short"] == jc["call_long"])), 3) if len(jc) else np.nan,
            "kappa_confident": round(cohen_kappa(jc["call_short"], jc["call_long"]), 3),
            "short_frac_tet": round(cs["tet"], 3),
            "long_frac_tet": round(cl["tet"], 3),
            "tet_shift_long_minus_short": round(cl["tet"] - cs["tet"], 3),
        })
    return pd.DataFrame(rows)


def plot_shift(pg, conc, outpath):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return
    if conc.empty:
        return
    colors = {"dip": "#4C72B0", "tri": "#DD8452", "tet": "#55A868"}
    labels, bars = [], []
    pgi = pg.set_index("sample")
    for _, r in conc.iterrows():
        suffix = "_ont" if r["long_platform"] == "ONT" else "_pb"
        for samp, lab in [(r["strain"], f"{r['strain']}\nshort"), (r["strain"] + suffix, f"{r['strain']}\n{r['long_platform']}")]:
            labels.append(lab)
            bars.append([pgi.loc[samp, f"frac_{p}"] for p in PLOIDIES])
    bars = np.array(bars, dtype=float)
    fig, ax = plt.subplots(figsize=(max(8, 0.8 * len(labels)), 5))
    x = np.arange(len(labels))
    bottom = np.zeros(len(labels))
    for i, p in enumerate(PLOIDIES):
        ax.bar(x, bars[:, i], bottom=bottom, color=colors[p], label={"dip": "2n", "tri": "3n", "tet": "4n"}[p])
        bottom += np.nan_to_num(bars[:, i])
    for k in range(1, len(labels) // 2):
        ax.axvline(2 * k - 0.5, color="grey", lw=0.6, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_ylabel("Fraction of called genes")
    ax.set_title("nQuire per-gene ploidy calls: short vs long reads, same strain")
    kap = "  ".join(f"{r.strain}: κ={r.kappa:.2f}" for r in conc.itertuples())
    fig.text(0.01, 0.005, kap, fontsize=7)
    ax.legend(title="winning model", bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(outpath)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", default="results/nquire2")
    ap.add_argument("--outdir", default=None, help="default: <indir>/concordance")
    ap.add_argument("--min-margin", type=float, default=5.0,
                    help="delta-LL margin for a 'confident' call (default 5)")
    a = ap.parse_args()
    outdir = a.outdir or os.path.join(a.indir, "concordance")
    os.makedirs(outdir, exist_ok=True)

    res = load_pergene(a.indir)
    pg = summarise_pergene(res, a.min_margin)
    pg.to_csv(os.path.join(outdir, "pergene_call_summary.tsv"), sep="\t", index=False)

    wg = load_wholegenome(a.indir)
    wg.to_csv(os.path.join(outdir, "wholegenome_summary.tsv"), sep="\t", index=False)

    conc = concordance(res, a.min_margin, outdir)
    conc.to_csv(os.path.join(outdir, "platform_concordance.tsv"), sep="\t", index=False)

    plot_shift(pg, conc, os.path.join(outdir, "platform_shift.pdf"))

    pd.set_option("display.width", 200)
    print(pg.to_string(index=False), "\n")
    print(wg[[c for c in ["sample", "platform", "status", "call", "margin", "relfit"] if c in wg]].to_string(index=False), "\n")
    print(conc.to_string(index=False))
    print(f"\nWrote outputs to {outdir}/")


if __name__ == "__main__":
    main()
