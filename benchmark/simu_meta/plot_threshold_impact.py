#!/usr/bin/env python3
"""
plot_threshold_impact.py — how each threshold's VALUE changes accuracy (Part G / C5).

Computed on the bg_300 community (300 genomes — the microbiome-scale complexity point):
x-axis = the parameter value, y-axis = recall / precision / F1 of the ECE-host assignment.
One parameter is swept while the others are held at default (score>0.5, specificity<0.01).
Denominator for recall = every planted ECE (unassigned counts as a miss).

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

SETTING = "bg_300"                       # microbiome-scale complexity point
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


def gather(label):
    """best-host calls (host_summary) + total planted-ECE count for one community."""
    d = f"{ROOT}/{label}"
    hs = f"{d}/modifi/{label}/host_summary.csv"
    mge = f"{d}/{label}.mge_list.tsv"
    if not (os.path.exists(hs) and os.path.exists(mge)):
        return None
    total = pd.read_csv(mge, sep="\t")["seq_name"].nunique()
    h = pd.read_csv(hs)
    h["correct"] = h["MGE"].map(sra) == h["host"].map(sra)
    h["ece_frac"] = h["motif_info"].apply(mean_ece_frac)
    return h[["MGE", "final_score", "specificity", "total_sites", "ece_frac", "correct"]], total


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
    g = gather(SETTING)
    if g is None:
        raise SystemExit(f"{SETTING} host_summary/mge_list not found")
    best, total = g
    print(f"[plot] {SETTING}: {len(best)} scored ECEs / {total} planted")

    fig, ax = plt.subplots(2, 2, figsize=(13, 10.5))
    for a, (param, grid, name, cited, xscale) in zip(ax.ravel(), PANELS):
        R, P, F = metrics_vs(best, total, param, grid)
        a.plot(grid, R, color=BLUE, lw=2, label="recall")
        a.plot(grid, P, color=ORANGE, lw=2, label="precision")
        a.plot(grid, F, color=GREEN, lw=2, label="F1")
        a.axvline({"final_score": S0, "specificity": P0, "ece_frac": 0.0, "total_sites": 0}[param],
                  color=GREY, ls=":", lw=1.2)
        if xscale == "log":
            a.set_xscale("log")
        a.set(xlabel=f"{name} threshold", ylabel="metric", ylim=(0, 1.10),
              title=f"Accuracy vs {name}")
        a.legend(loc="lower left", fontsize=9)
    fig.suptitle(f"How each threshold changes ECE-host accuracy — {SETTING} "
                 f"(300 genomes, microbiome-scale complexity; {total} planted ECEs)",
                 fontsize=12.5, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    pdf = f"{OUT}/threshold_impact.pdf"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(pdf.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"[plot] wrote {pdf}")


if __name__ == "__main__":
    main()
