#!/usr/bin/env python3
"""plot_e2e_completeness.py — ECE->host linkage success vs host-genome assembly completeness.

Pools the de-novo e2e rep1 communities (bg_80/150/300 + toy2) where each ECE's true host
genome has a CheckM2 completeness. Bins by completeness and plots the fraction of ECEs
correctly linked (at the default operating point). Shows that incompletely assembled hosts
lose their ECE links, while completely assembled hosts are reliably recovered.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/mag_completeness"
COMMS = ["bg_80", "bg_150", "bg_300"]
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
BLUE = "#0072B2"


def main():
    frames = []
    for c in COMMS:
        p = f"{ROOT}/C4_{c}/e2e_predictions.csv"
        if not os.path.exists(p):
            continue
        d = pd.read_csv(p)
        d["ok"] = ((d.final_score > 0.5) & (d.specificity < 0.01) &
                   (d.ece_origin == d.host_origin) & d.ece_origin.notna() &
                   (d.ece_origin != "UNMAPPED"))
        frames.append(d[["host_completeness", "ok", "final_score"]])
    a = pd.concat(frames, ignore_index=True).dropna(subset=["host_completeness"])

    bins = [0, 60, 70, 80, 90, 95, 100.01]
    labs = ["<60", "60-70", "70-80", "80-90", "90-95", "95-100"]
    a["bin"] = pd.cut(a.host_completeness, bins, labels=labs)
    g = a.groupby("bin", observed=True).agg(n=("ok", "size"), correct=("ok", "sum"))
    g["success"] = g.correct / g.n
    # Wilson CI
    z = 1.96
    p, n = g.success.values, g.n.values
    denom = 1 + z**2 / n
    ci = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    g["ci"] = ci

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    # A: success rate by completeness bin
    xs = range(len(g))
    ax[0].bar(xs, g.success, yerr=g.ci, color=BLUE, width=0.7, capsize=3)
    for i, (s, nn) in enumerate(zip(g.success, g.n)):
        ax[0].annotate(f"{s:.2f}\n(n={nn})", (i, s), textcoords="offset points",
                       xytext=(0, 4), ha="center", fontsize=8.5)
    ax[0].set_xticks(list(xs)); ax[0].set_xticklabels(g.index.astype(str))
    ax[0].set(ylim=(0, 1.1), xlabel="host-genome CheckM2 completeness (%)",
              ylabel="fraction of ECEs correctly linked",
              title="A. Linkage success vs host completeness")
    # B: per-ECE scatter (completeness vs linkage score, colored by correct)
    col = a.ok.map({True: "#0072B2", False: "#D55E00"})
    ax[1].scatter(a.host_completeness, a.final_score, c=col, s=22, edgecolor="#555555", lw=0.3, alpha=0.8)
    ax[1].axhline(0.5, ls="--", color="#888888", lw=1)
    ax[1].set(xlabel="host-genome completeness (%)", ylabel="ECE linkage score",
              title="B. Per-ECE (blue=correct, orange=miss/declined)")
    fig.suptitle("ECE->host linkage vs assembly completeness (de-novo, pooled bg_80/150/300; rep1)",
                 fontsize=12.5, y=1.02)
    fig.tight_layout()
    out = f"{OUT}/e2e_linkage_vs_completeness.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    print(g.to_string())


if __name__ == "__main__":
    main()
