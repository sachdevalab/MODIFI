#!/usr/bin/env python3
"""plot_orphan_calibration.py — threshold calibration WITH a real negative class (orphans).

Uses bg80_orph, where 37 ECEs are orphans (true host removed) and 32 are planted. Sweeps
the two linkage-decision thresholds (specificity, final_score) and plots recall, precision,
and orphan-FPR — so the precision/FPR tradeoff is grounded in host-absent negatives, not
asserted. Vector PDF -> tmp/rev_figs/simu_meta/threshold/orphan_calibration.pdf
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/home/shuaiw/borg/paper/simu_meta_dir/C1/bg80_orph"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE, RED, GREY = "#0072B2", "#D55E00", "#CC3311", "#888888"
S0, P0 = 0.5, 0.01


def sra(x):
    return str(x).split("_")[0]


def load():
    orph = set(l.strip() for l in open(f"{D}/bg80_orph.orphans.txt"))
    h = pd.read_csv(f"{D}/modifi/bg80_orph/host_summary.csv")
    h["is_orphan"] = h["MGE"].isin(orph)
    h["correct"] = (~h.is_orphan) & (h["MGE"].map(sra) == h["host"].map(sra))
    return h, int((~h.is_orphan).sum()), int(h.is_orphan.sum())


def metrics(h, n_pl, n_or, s, p):
    acc = h[(h.final_score > s) & (h.specificity < p)]
    tp = int(((~acc.is_orphan) & acc.correct).sum())
    orphFP = int(acc.is_orphan.sum())
    n = len(acc)
    return tp / n_pl, (tp / n if n else np.nan), orphFP / n_or


def main():
    h, n_pl, n_or = load()
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

    # panel A: sweep specificity (final_score held at default)
    pgrid = np.logspace(-4, 0, 80)
    R, P, F = zip(*[metrics(h, n_pl, n_or, S0, p) for p in pgrid])
    ax[0].plot(pgrid, R, color=BLUE, lw=2, label="recall")
    ax[0].plot(pgrid, P, color=ORANGE, lw=2, label="precision")
    ax[0].plot(pgrid, F, color=RED, lw=2, ls="--", label="orphan-FPR")
    for xv, lab in [(0.01, "default\n0.01"), (0.001, "cited\n0.001")]:
        ax[0].axvline(xv, color=GREY, ls=":", lw=1)
        ax[0].text(xv, 1.02, lab, ha="center", va="bottom", fontsize=7.5, color=GREY)
    ax[0].set_xscale("log")
    ax[0].set(xlabel="specificity accept threshold (< x)", ylabel="metric", ylim=(0, 1.10),
              title=f"Accuracy vs specificity  (final_score>{S0})")
    ax[0].legend(loc="center left", fontsize=9)

    # panel B: sweep final_score (specificity held at default)
    sgrid = np.linspace(0, 1, 80)
    R, P, F = zip(*[metrics(h, n_pl, n_or, s, P0) for s in sgrid])
    ax[1].plot(sgrid, R, color=BLUE, lw=2, label="recall")
    ax[1].plot(sgrid, P, color=ORANGE, lw=2, label="precision")
    ax[1].plot(sgrid, F, color=RED, lw=2, ls="--", label="orphan-FPR")
    for xv, lab in [(0.5, "default\n0.5"), (0.8, "cited\n0.8")]:
        ax[1].axvline(xv, color=GREY, ls=":", lw=1)
        ax[1].text(xv, 1.02, lab, ha="center", va="bottom", fontsize=7.5, color=GREY)
    ax[1].set(xlabel="final_score accept threshold (> x)", ylabel="metric", ylim=(0, 1.10),
              title=f"Accuracy vs linkage score  (specificity<{P0})")
    ax[1].legend(loc="center left", fontsize=9)

    fig.suptitle(f"Threshold calibration with host-absent orphans as negatives "
                 f"(bg80_orph: {n_pl} planted + {n_or} orphans)", fontsize=12, y=1.02)
    fig.tight_layout()
    out = f"{OUT}/orphan_calibration.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
