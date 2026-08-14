#!/usr/bin/env python3
"""score_mag.py — score the synthetic completeness x contamination grid on fragmented ladder_58.

Each cell has a host_summary (ECE -> best host bin, produced by --run_steps host with the
cell's --bin_file). An ECE link is CORRECT iff its best host bin (an isolate) == the ECE's
true isolate (SRA prefix of the ECE contig). Reports recall/precision per (completeness,
contamination) and writes a heatmap + the per-cell table.

Usage: score_mag.py <results_dir>   (dir of c<comp>_x<contam>.host_summary.csv files)
"""
import os, sys, re, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S0, P0 = 0.5, 0.01


def sra(x):
    return str(x).split("_")[0] if pd.notna(x) else None


def main():
    RES = sys.argv[1]
    rows = []
    for f in sorted(glob.glob(os.path.join(RES, "c*_x*.host_summary.csv"))):
        m = re.search(r"c(\d+)_x(\d+)\.host_summary", os.path.basename(f))
        if not m:
            continue
        comp, contam = int(m.group(1)), int(m.group(2))
        h = pd.read_csv(f)
        # total planted ECEs = the mge denominator (host_summary rows = ECEs with a best host)
        h["assigned"] = (h["final_score"] > S0) & (h["specificity"] < P0)
        h["correct"] = h["assigned"] & (h["MGE"].map(sra) == h["host"].map(sra))
        total = len(h)
        tp = int(h["correct"].sum())
        n_assigned = int(h["assigned"].sum())
        rows.append(dict(completeness=comp, contamination=contam, n_ece=total,
                         tp=tp, assigned=n_assigned,
                         recall=tp / total if total else np.nan,
                         precision=tp / n_assigned if n_assigned else np.nan))
    d = pd.DataFrame(rows).sort_values(["contamination", "completeness"])
    d.to_csv(os.path.join(RES, "frag_grid_summary.csv"), index=False)
    print(d.to_string(index=False))

    # heatmap of recall over the completeness x contamination grid
    if len(d):
        piv = d.pivot(index="contamination", columns="completeness", values="recall")
        fig, ax = plt.subplots(figsize=(7, 4.2))
        im = ax.imshow(piv.values, aspect="auto", cmap="viridis", vmin=0, vmax=1, origin="lower")
        ax.set_xticks(range(len(piv.columns))); ax.set_xticklabels(piv.columns)
        ax.set_yticks(range(len(piv.index))); ax.set_yticklabels(piv.index)
        ax.set(xlabel="MAG completeness (%)", ylabel="contamination (%)",
               title="ECE->host recall vs MAG completeness & contamination\n(fragmented ladder_58, per-MAG bins)")
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            color="white" if v < 0.6 else "black", fontsize=9)
        fig.colorbar(im, label="recall")
        fig.tight_layout()
        out = os.path.join(RES, "frag_completeness_contamination.pdf")
        fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150)
        print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
