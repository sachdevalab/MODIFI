#!/usr/bin/env python3
"""plot_c3_coverage.py — ECE-host recall/precision vs donor sequencing depth (C3 titration).

Fixed bg_80 composition (40 donor species + 40 background) re-run at donor depth
{5,10,20,30,40x} (background fixed 10x). Locates the usable-depth floor: recall collapses
at 5x (MODIFI's min_ctg_cov floor) and saturates by ~30x, while precision stays ~1.0.
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/coverage"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE = "#0072B2", "#D55E00"
POINTS = [(10, "cov_d10"), (20, "cov_d20"), (30, "bg_80"), (40, "cov_d40")]


def sra(x):
    return str(x).split("_")[0]


def main():
    dep, rec, prec = [], [], []
    for d, lbl in POINTS:
        f = f"{ROOT}/{lbl}/modifi/{lbl}/host_summary.csv"
        if not os.path.exists(f):
            continue
        h = pd.read_csv(f)
        tot = pd.read_csv(f"{ROOT}/{lbl}/{lbl}.mge_list.tsv", sep="\t").seq_name.nunique()
        h["p"] = (h.final_score > 0.5) & (h.specificity < 0.01)
        tp = (h["p"] & (h.MGE.map(sra) == h.host.map(sra))).sum()
        na = h["p"].sum()
        dep.append(d); rec.append(tp / tot); prec.append(tp / na if na else 0.0)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(dep, rec, "-o", color=BLUE, lw=2, ms=8, label="recall")
    ax.plot(dep, prec, "-s", color=ORANGE, lw=2, ms=7, label="precision")
    for x, r, p in zip(dep, rec, prec):
        ax.annotate(f"{r:.2f}", (x, r), textcoords="offset points", xytext=(0, -15), ha="center", fontsize=9, color=BLUE)
        ax.annotate(f"{p:.2f}", (x, p), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=9, color=ORANGE)
    ax.set(xlabel="donor sequencing depth (x)", ylabel="metric", ylim=(0, 1.05),
           title="C3: ECE-host recall & precision vs coverage\n(bg_80 composition; background fixed 10x; rep1)")
    ax.set_xticks(dep)
    ax.legend(loc="center right", fontsize=10)
    fig.tight_layout()
    out = f"{OUT}/c3_coverage_titration.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out} | depth->recall: {dict(zip(dep, [round(r,2) for r in rec]))}")


if __name__ == "__main__":
    main()
