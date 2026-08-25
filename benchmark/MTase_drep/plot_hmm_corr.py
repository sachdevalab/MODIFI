#!/usr/bin/env python
"""
Original HMM-group MTase count vs No. of motifs correlation (kept as the baseline
companion to the mmseqs-dereplicated figure). MT_hmm_num is identical across the
cutoff sweep, so we read it from one sourcedata file.

Outputs (tmp/rev_figs/MTase_drep/):
  MTase_hmm_vs_motif_num.{pdf,png}
  mtase_hmm_motif_corr.sourcedata.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, linregress

OUT_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/MTase_drep"
SRC = os.path.join(OUT_DIR, "mtase_drep_motif_corr.id90_cov80.sourcedata.csv")


def main():
    df = pd.read_csv(SRC)[
        ["sample", "environment", "contig", "motif_num", "MT_hmm_num", "ctg_len"]].copy()
    df.to_csv(os.path.join(OUT_DIR, "mtase_hmm_motif_corr.sourcedata.csv"), index=False)

    r, p = pearsonr(df["MT_hmm_num"], df["motif_num"])
    lr = linregress(df["MT_hmm_num"], df["motif_num"])

    envs = sorted(df["environment"].dropna().unique())
    palette = dict(zip(envs, sns.color_palette("tab10", len(envs))))

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    rng = np.random.default_rng(0)
    jx = df["MT_hmm_num"] + rng.uniform(-0.18, 0.18, len(df))
    jy = df["motif_num"] + rng.uniform(-0.18, 0.18, len(df))
    for env in envs:
        m = df["environment"] == env
        ax.scatter(jx[m], jy[m], s=16, alpha=0.65, color=palette[env],
                   edgecolor="none", label=env)
    sns.regplot(x="MT_hmm_num", y="motif_num", data=df, scatter=False,
                ax=ax, color="black", line_kws={"linewidth": 1.4}, truncate=False)
    ax.set_xlabel("No. of MTases (HMM group)")
    ax.set_ylabel("No. of motifs")
    ax.set_title(f"Pearson r = {r:.2f}, p = {p:.2e}  (n = {len(df)})")
    ax.legend(title="Environment", fontsize=7, title_fontsize=8,
              markerscale=1.4, frameon=False, loc="best")
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "MTase_hmm_vs_motif_num.pdf"))
    fig.savefig(os.path.join(OUT_DIR, "MTase_hmm_vs_motif_num.png"), dpi=300)

    print("=== Original HMM-group MTase vs No. of motifs (n = %d) ===" % len(df))
    print(f"  Pearson r = {r:.4f}, p = {p:.3e}")
    print(f"  slope = {lr.slope:.4f} (SE {lr.stderr:.4f}), intercept = {lr.intercept:.4f}")
    print(f"  mean HMM-group MTase per genome = {df['MT_hmm_num'].mean():.3f}")
    print(f"\nWrote:\n  {os.path.join(OUT_DIR, 'MTase_hmm_vs_motif_num.pdf')}"
          f"\n  {os.path.join(OUT_DIR, 'mtase_hmm_motif_corr.sourcedata.csv')}")


if __name__ == "__main__":
    main()
