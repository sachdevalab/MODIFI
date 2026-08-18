#!/usr/bin/env python3
"""
ECE-host linkage across a gradient of host MAG/contig completeness.

Reviewer response: assess ECE-host linkage performance across a gradient of MAG
completeness. For each of the real run2 samples, we take the contigs MODIFI
actually profiled for methylation (the candidate host universe), join their
CheckM2 completeness, and mark which ones received a confident ECE linkage
(default filter: specificity < 0.01 AND final_score > 0.5, as in
scripts/estimate_linkage.py). Pooled across samples we plot:
  A. linkage rate  = fraction of candidate host contigs linked to >=1 ECE, per
                     completeness bin (Wilson 95% CI)
  B. linkage count = number of ECE-host linkages, per completeness bin

A secondary variant restricts to contamination <= 5% (the paper criterion).
"""

import os
import glob
import math

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})

# Okabe-Ito colorblind-safe palette (repo house style)
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#888888"

RUN2 = "/home/shuaiw/borg/paper/run2"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/complete_assess"

# MODIFI default linkage filter (scripts/estimate_linkage.py:587-588)
SPEC_CUT = 0.01      # specificity < 0.01
SCORE_CUT = 0.5      # final_score > 0.5

# Completeness bins (contamination ignored for this axis).
COMP_EDGES = [0, 25, 50, 75, 90, 100.01]
COMP_LABELS = ["0-25", "25-50", "50-75", "75-90", "90-100"]

# Contamination bins (completeness ignored for this axis).
CONT_EDGES = [0, 1, 2, 5, 10, 1e9]
CONT_LABELS = ["0-1", "1-2", "2-5", "5-10", ">10"]


def wilson_ci(k, n, z=1.96):
    """Wilson score 95% CI for a binomial proportion; returns (lo, hi)."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def build_contig_table():
    """Return a pooled per-contig DataFrame over all profiled contigs across
    samples, with completeness, contamination, linked flag and n_linkages."""
    rows = []
    n_samples = 0
    n_missing_checkm = 0
    meth_dirs = sorted(glob.glob(os.path.join(RUN2, "*", "*_methylation4")))
    for mdir in meth_dirs:
        sample_dir = os.path.dirname(mdir)
        sample = os.path.basename(sample_dir)

        contigs_list = os.path.join(mdir, "contigs_list.txt")
        host_summary = os.path.join(mdir, "host_summary.csv")
        checkm = os.path.join(sample_dir, "checkM2", "quality_report.tsv")
        if not (os.path.exists(contigs_list) and os.path.exists(checkm)):
            continue
        n_samples += 1

        # Candidate host universe: contigs MODIFI profiled for methylation.
        profiled = []
        with open(contigs_list) as fh:
            for line in fh:
                p = line.strip()
                if not p:
                    continue
                cid = os.path.basename(p)
                if cid.endswith(".fa"):
                    cid = cid[:-3]
                profiled.append(cid)
        profiled = set(profiled)

        # CheckM2 completeness / contamination per contig.
        cm = pd.read_csv(checkm, sep="\t")
        cm = cm[["Name", "Completeness", "Contamination"]].copy()
        comp = dict(zip(cm["Name"], cm["Completeness"]))
        cont = dict(zip(cm["Name"], cm["Contamination"]))

        # Confident linkages (default filter). host_summary may be empty/absent.
        linked_counts = {}
        if os.path.exists(host_summary):
            hs = pd.read_csv(host_summary)
            if len(hs) > 0 and {"specificity", "final_score", "host"} <= set(hs.columns):
                hs = hs[(hs["specificity"] < SPEC_CUT) & (hs["final_score"] > SCORE_CUT)]
                linked_counts = hs["host"].value_counts().to_dict()

        for cid in profiled:
            if cid not in comp:
                n_missing_checkm += 1
                continue
            nlink = int(linked_counts.get(cid, 0))
            rows.append({
                "sample": sample,
                "contig": cid,
                "completeness": float(comp[cid]),
                "contamination": float(cont[cid]),
                "linked": nlink > 0,
                "n_linkages": nlink,
            })

    df = pd.DataFrame(rows)
    print(f"pooled {n_samples} samples, {len(df)} profiled contigs with CheckM2 "
          f"({n_missing_checkm} profiled contigs dropped for missing CheckM2 row)")
    return df


def summarize(df, col, edges, labels):
    """Per-bin summary DataFrame, binning contigs on `col`."""
    df = df.copy()
    df["bin"] = pd.cut(df[col], bins=edges, labels=labels,
                       right=False, include_lowest=True)
    out = []
    for i, lab in enumerate(labels):
        sub = df[df["bin"] == lab]
        n_contigs = len(sub)
        n_linked = int(sub["linked"].sum())
        n_linkages = int(sub["n_linkages"].sum())
        rate = n_linked / n_contigs if n_contigs else 0.0
        lo, hi = wilson_ci(n_linked, n_contigs)
        out.append({
            "bin": lab,
            "low": edges[i],
            "high": edges[i + 1],
            "n_contigs": n_contigs,
            "n_linked": n_linked,
            "linkage_rate": rate,
            "ci_low": lo,
            "ci_high": hi,
            "n_linkages": n_linkages,
        })
    return pd.DataFrame(out)


def rate_panel(ax, summary, labels, xlabel, title):
    """Draw a linkage-rate barplot (Wilson CI) for one binning on `ax`."""
    xs = np.arange(len(labels))
    rate = summary["linkage_rate"].values
    lo = summary["ci_low"].values
    hi = summary["ci_high"].values
    yerr = np.vstack([rate - lo, hi - rate])
    ncontig = summary["n_contigs"].values

    ax.bar(xs, rate, yerr=yerr, color=BLUE, width=0.7, capsize=3,
           edgecolor=GREY, linewidth=0.5)
    for x, r, n in zip(xs, rate, ncontig):
        ax.annotate(f"{r:.3f}\n(n={n})", (x, hi[x]), textcoords="offset points",
                    xytext=(0, 4), ha="center", va="bottom", fontsize=8.5)
    ax.set(xticks=xs, ylim=(0, 1.18), xlabel=xlabel,
           ylabel="Fraction of host contigs linked to ≥1 ECE", title=title)
    ax.set_xticklabels(labels)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = build_contig_table()
    n_samp = df["sample"].nunique()

    # Bin by completeness (contamination ignored) and by contamination
    # (completeness ignored).
    summ = summarize(df, "completeness", COMP_EDGES, COMP_LABELS)
    summ.to_csv(os.path.join(OUT, "linkage_vs_completeness_summary.tsv"),
                sep="\t", index=False)
    print(summ.to_string(index=False))

    # Contamination axis: restrict to reasonably complete contigs (completeness >= 50%)
    # so the contamination effect is not confounded by incompleteness. Without this, the
    # bins are dominated by the many very-incomplete contigs and contamination is
    # meaningless (a contig cannot be assessed for contamination when barely assembled).
    COMP_FLOOR = 50.0
    df_c = df[df["completeness"] >= COMP_FLOOR]
    summ_c = summarize(df_c, "contamination", CONT_EDGES, CONT_LABELS)
    summ_c.to_csv(os.path.join(OUT, "linkage_vs_contamination_summary.tsv"),
                  sep="\t", index=False)
    print(f"\n[by contamination, completeness >= {COMP_FLOOR:.0f}%: "
          f"{len(df_c)}/{len(df)} contigs]")
    print(summ_c.to_string(index=False))

    # Combined figure: A = completeness, B = contamination (both linkage rate).
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    rate_panel(ax[0], summ, COMP_LABELS, "Host contig completeness (%)",
               "a. Linkage rate vs. host completeness")
    rate_panel(ax[1], summ_c, CONT_LABELS, "Host contig contamination (%)",
               f"b. Linkage rate vs. host contamination\n(completeness ≥ {COMP_FLOOR:.0f}%)")
    fig.tight_layout()
    out_pdf = os.path.join(OUT, "linkage_vs_mag_quality.pdf")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_pdf.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
