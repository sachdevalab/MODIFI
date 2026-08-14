#!/usr/bin/env python3
"""plot_c1_complexity.py — species-level recall/precision vs community complexity (C1).

Donor-only ladder (25/40/58 species) + background-scaled communities (80/150/300 genomes).
Mean +/- 95% CI over replicates. Shows complexity's impact: precision ~invariant, recall
flat to ~150 genomes then a modest dip at 300 (borderline ECEs conservatively declined).
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AGG = None
for root, _, files in os.walk("/home/shuaiw/borg/paper/simu_meta_dir/C1"):
    if "size_summary_agg.csv" in files and "/eval" in root:
        AGG = os.path.join(root, "size_summary_agg.csv"); break
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/complexity"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
KEEP = ["ladder_25", "ladder_40", "ladder_58", "bg_80", "bg_150", "bg_300"]
BLUE, ORANGE = "#0072B2", "#D55E00"


def main():
    d = pd.read_csv(AGG)
    d = d[d.setting.isin(KEEP)].copy().sort_values("size")
    fig, ax = plt.subplots(figsize=(8, 5))
    for col, c, lab in [("recall", BLUE, "recall"), ("precision", ORANGE, "precision")]:
        m = d[f"{col}_mean"]; ci = d[f"{col}_ci"].fillna(0)
        ax.plot(d["size"], m, "-o", color=c, lw=2, ms=7, label=lab)
        ax.fill_between(d["size"], (m - ci).clip(0, 1), (m + ci).clip(0, 1), color=c, alpha=0.18, lw=0)
    # mark donor-only vs background regions
    ax.axvspan(20, 58, color="#EEEEEE", alpha=0.5, zorder=0)
    ax.text(40, 0.02, "donor-only ladder", ha="center", fontsize=8, color="#666")
    ax.text(190, 0.02, "+ background (microbiome scale)", ha="center", fontsize=8, color="#666")
    for _, r in d.iterrows():
        ax.annotate(f"{r['recall_mean']:.2f}", (r["size"], r["recall_mean"]),
                    textcoords="offset points", xytext=(0, 9), ha="center", fontsize=8, color=BLUE)
    ax.set(xlabel="community size (genomes)", ylabel="species-level metric",
           ylim=(0, 1.05), title="Complexity impact on MODIFI ECE->host linkage (species level)")
    ax.set_xscale("log"); ax.set_xticks(d["size"]); ax.set_xticklabels(d["size"].astype(int))
    ax.legend(loc="lower left", fontsize=10)
    fig.tight_layout()
    out = f"{OUT}/c1_complexity_species.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out}  (n_rep: {dict(zip(d.setting, d.n_rep))})")


if __name__ == "__main__":
    main()
