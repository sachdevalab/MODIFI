#!/usr/bin/env python3
"""
plot_threshold_impact.py — how each threshold's VALUE changes accuracy (Part G / C5).

Directly interpretable companion to the ROC figure: x-axis = the parameter value,
y-axis = recall / precision / F1 of the ECE-host assignment, pooled over ALL completed
ladder + background communities. One parameter is swept while the others are held at
default (score>0.5, specificity<0.01). Replicates
(seeds 42-46) give a mean +/- SD shaded band. Denominator for recall = every planted
ECE (unassigned counts as a miss).

Vector PDF -> /home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold/threshold_impact.pdf
"""
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#888888"

SETTINGS = ["ladder_25", "ladder_40", "ladder_58", "bg_80", "bg_150", "bg_300"]
S0, P0 = 0.5, 0.01                       # default holds while sweeping another parameter


def sra(x): return str(x).split("_")[0]


def mean_ece_frac(mi):
    if not isinstance(mi, str) or not mi:
        return 0.0
    fr = []
    for e in mi.split(";"):
        p = e.split(":")
        if len(p) >= 6:
            try:
                pt, pm = int(p[4]), int(p[5])
                if pt > 0:
                    fr.append(pm / pt)
            except ValueError:
                pass
    return float(np.mean(fr)) if fr else 0.0


def gather_rep(rep):
    """Pool best-host calls (host_summary) + total ECE count over all settings for one
    replicate. rep 1 = bare labels; rep>=2 = _repN."""
    rows, total = [], 0
    for s in SETTINGS:
        label = s if rep == 1 else f"{s}_rep{rep}"
        d = f"{ROOT}/{label}"
        hs = f"{d}/modifi/{label}/host_summary.csv"
        mge = f"{d}/{label}.mge_list.tsv"
        if not (os.path.exists(hs) and os.path.exists(mge)):
            continue
        total += len(pd.read_csv(mge, sep="\t"))
        h = pd.read_csv(hs)
        if h.empty:
            continue
        h["correct"] = h["MGE"].map(sra) == h["host"].map(sra)
        h["ece_frac"] = h["motif_info"].apply(mean_ece_frac)
        rows.append(h[["MGE", "final_score", "specificity", "total_sites", "ece_frac", "correct"]])
    if not rows or total == 0:
        return None
    return pd.concat(rows, ignore_index=True), total


def metrics_vs(best, total, param, grid):
    """recall/precision/F1 as `param`'s accept-threshold sweeps `grid`; others at default."""
    rec, prec, f1 = [], [], []
    for t in grid:
        s, p, k, f = S0, P0, 0.0, 0.0
        if param == "final_score": s = t
        elif param == "specificity": p = t
        elif param == "total_sites": k = t
        elif param == "ece_frac": f = t
        acc = best[(best.final_score > s) & (best.specificity < p) &
                   (best.total_sites >= k) & (best.ece_frac > f)]
        tp = int(acc.correct.sum()); n = len(acc)
        r = tp / total
        pr = tp / n if n else np.nan
        rec.append(r); prec.append(pr)
        f1.append(2 * pr * r / (pr + r) if (pr and r) else 0.0)
    return np.array(rec), np.array(prec), np.array(f1)


def band(ax, grid, mat, color, label):
    """mat: (n_rep, len(grid)); plot mean line + SD band."""
    m = np.nanmean(mat, 0)
    ax.plot(grid, m, color=color, lw=2, label=label)
    if mat.shape[0] > 1:
        sd = np.nanstd(mat, 0)
        ax.fill_between(grid, np.clip(m - sd, 0, 1), np.clip(m + sd, 0, 1), color=color, alpha=0.18, lw=0)


PANELS = [
    ("final_score", np.round(np.linspace(0, 1, 51), 3), "linkage score", 0.8, "linear"),
    ("specificity", np.logspace(-4, 0, 51), "specificity", 0.001, "log"),
    ("ece_frac", np.round(np.linspace(0, 1, 51), 3), "motif methylation fraction", 0.3, "linear"),
    ("total_sites", np.round(np.linspace(0, 400, 51)), "# modified sites", 100, "linear"),
]


def main():
    reps = [r for r in range(1, 6) if gather_rep(r)]
    data = {r: gather_rep(r) for r in reps}
    print(f"[plot] replicates with data: {reps}")

    fig, ax = plt.subplots(2, 2, figsize=(13, 10.5))
    for a, (param, grid, name, cited, xscale) in zip(ax.ravel(), PANELS):
        R = np.vstack([metrics_vs(*data[r], param, grid)[0] for r in reps])
        P = np.vstack([metrics_vs(*data[r], param, grid)[1] for r in reps])
        F = np.vstack([metrics_vs(*data[r], param, grid)[2] for r in reps])
        band(a, grid, R, BLUE, "recall")
        band(a, grid, P, ORANGE, "precision")
        band(a, grid, F, GREEN, "F1")
        if xscale == "log":
            a.set_xscale("log")
        a.set(xlabel=f"{name} threshold", ylabel="metric", ylim=(0, 1.10),
              title=f"Accuracy vs {name}")
        a.legend(loc="lower left", fontsize=9)
    n = len(reps)
    fig.suptitle(f"How each threshold changes ECE-host accuracy "
                 f"(pooled ladder+background; band = ±SD over {n} replicate"
                 f"{'s' if n>1 else ''})", fontsize=12.5, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    pdf = f"{OUT}/threshold_impact.pdf"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(pdf.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"[plot] wrote {pdf}")


if __name__ == "__main__":
    main()
