#!/usr/bin/env python
"""
Aggregate the MTase-dereplication cutoff sweep and help pick the best threshold.

Reads every mtase_drep_motif_corr.id*_cov*.sourcedata.csv produced by
mtase_drep_corr.py across identity cutoffs, then writes:
  - cutoff_comparison.sourcedata.csv : one row per cutoff (Pearson r/p vs motif_num,
    Spearman rho, mean/median unique MTase count, total clusters, ...)
  - MTase_drep_cutoff_comparison.{pdf,png} : (left) Pearson r vs identity cutoff,
    (right) small-multiple scatter of unique-MTase vs motif_num per cutoff.

Run after the SLURM array (run_mtase_drep_slurm.sh) finishes.
"""

import os
import glob
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr, linregress

OUT_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/MTase_drep"


def parse_tag(path):
    m = re.search(r"\.id(\d+)_cov(\d+)\.sourcedata\.csv$", os.path.basename(path))
    return (int(m.group(1)) / 100.0, int(m.group(2)) / 100.0) if m else (None, None)


def main():
    files = sorted(glob.glob(os.path.join(
        OUT_DIR, "mtase_drep_motif_corr.id*_cov*.sourcedata.csv")))
    if not files:
        raise SystemExit("No per-cutoff sourcedata found; run the SLURM array first.")

    rows = []
    frames = []  # (id, cov, df) for the small-multiple panel
    for f in files:
        sid, scov = parse_tag(f)
        df = pd.read_csv(f)
        r, p = pearsonr(df["MT_drep_num"], df["motif_num"])
        rho, prho = spearmanr(df["MT_drep_num"], df["motif_num"])
        r_hmm, p_hmm = pearsonr(df["MT_hmm_num"], df["motif_num"])
        # OLS regression of No. of motifs (y) on unique MTase count (x)
        lr = linregress(df["MT_drep_num"], df["motif_num"])
        lr_hmm = linregress(df["MT_hmm_num"], df["motif_num"])
        rows.append({
            "min_seq_id": sid, "min_cov": scov, "n_contigs": len(df),
            "pearson_r": r, "pearson_p": p, "spearman_rho": rho, "spearman_p": prho,
            "slope": lr.slope, "intercept": lr.intercept, "slope_stderr": lr.stderr,
            "pearson_r_hmm": r_hmm, "slope_hmm": lr_hmm.slope,
            "intercept_hmm": lr_hmm.intercept,
            "mean_unique_mtase": df["MT_drep_num"].mean(),
            "median_unique_mtase": df["MT_drep_num"].median(),
            "max_unique_mtase": int(df["MT_drep_num"].max()),
            "total_unique_mtase": int(df["MT_drep_num"].sum()),
            "mean_hmm_mtase": df["MT_hmm_num"].mean(),
        })
        frames.append((sid, scov, df))

    comp = pd.DataFrame(rows).sort_values("min_seq_id").reset_index(drop=True)
    comp_path = os.path.join(OUT_DIR, "cutoff_comparison.sourcedata.csv")
    comp.to_csv(comp_path, index=False)

    # ---- figure: r vs cutoff (left) + small-multiple scatters (right) ----
    frames.sort(key=lambda t: t[0])
    n = len(frames)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(4 + 2.4 * ncol, max(4.5, 2.4 * nrow)))
    gs = fig.add_gridspec(nrow, ncol + 2)

    axl = fig.add_subplot(gs[:, :2])
    axl.plot(comp["min_seq_id"], comp["pearson_r"], "-o", color="#1f77b4",
             label="mmseqs unique MTase")
    axl.axhline(comp["pearson_r_hmm"].iloc[0], color="grey", ls="--", lw=1,
                label="HMM-group count")
    axl.set_xlabel("mmseqs identity cutoff (coverage 0.8)")
    axl.set_ylabel("Pearson r (unique MTase vs No. of motifs)")
    axl.set_title("Correlation vs identity cutoff")
    axl.legend(fontsize=8, frameon=False)
    sns.despine(ax=axl)

    rng = np.random.default_rng(0)
    for i, (sid, scov, df) in enumerate(frames):
        ax = fig.add_subplot(gs[i // ncol, 2 + i % ncol])
        jx = df["MT_drep_num"] + rng.uniform(-0.18, 0.18, len(df))
        jy = df["motif_num"] + rng.uniform(-0.18, 0.18, len(df))
        ax.scatter(jx, jy, s=6, alpha=0.4, color="#444444", edgecolor="none")
        r, p = pearsonr(df["MT_drep_num"], df["motif_num"])
        ax.set_title(f">={int(sid*100)}% id  r={r:.2f}", fontsize=8)
        ax.tick_params(labelsize=6)
        if i // ncol == nrow - 1:
            ax.set_xlabel("unique MTase", fontsize=7)
        if i % ncol == 0:
            ax.set_ylabel("motifs", fontsize=7)
        sns.despine(ax=ax)

    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "MTase_drep_cutoff_comparison.pdf"))
    fig.savefig(os.path.join(OUT_DIR, "MTase_drep_cutoff_comparison.png"), dpi=300)

    print("=== MTase dereplication cutoff comparison (coverage 0.8) ===")
    print("Regression: No. of motifs = slope * (unique MTase) + intercept\n")
    focus = comp[["min_seq_id", "n_contigs", "pearson_r", "pearson_p",
                  "slope", "slope_stderr", "intercept",
                  "mean_unique_mtase"]].copy()
    with pd.option_context("display.float_format", lambda v: f"{v:.4g}"):
        print(focus.to_string(index=False))
    print(f"\nHMM-group baseline: r = {comp['pearson_r_hmm'].iloc[0]:.4f}, "
          f"slope = {comp['slope_hmm'].iloc[0]:.4f}, "
          f"intercept = {comp['intercept_hmm'].iloc[0]:.4f}")
    print("\nFull table:")
    print(comp.to_string(index=False))
    print(f"\nWrote:\n  {comp_path}"
          f"\n  {os.path.join(OUT_DIR, 'MTase_drep_cutoff_comparison.pdf')}")


if __name__ == "__main__":
    main()
