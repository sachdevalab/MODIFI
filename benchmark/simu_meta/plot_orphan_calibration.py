#!/usr/bin/env python3
"""plot_orphan_calibration.py — data-driven threshold calibration WITH a real negative
class (host-absent orphans), on orphan_300 x5 replicates.

Each community: 29 donor species keep their host (planted ECEs) and 29 have the host
removed (orphan ECEs, true host absent), over 242 background. This gives ground-truthed
positives AND negatives, so the two linkage-decision thresholds (final_score, specificity)
can be calibrated instead of asserted. Definitions at any operating point (s, p):
  pass       = final_score > s AND specificity < p
  recall     = planted-correct / n_planted
  precision  = planted-correct / all-confident (planted-correct + planted-mislink + orphan-FP)
  orphan-FPR = orphan-confident / n_orphan   (the host-absent false-positive rate)
One panel per hyperparameter (mean +/- 95% CI across 5 reps), the other three held at
default (final_score>0.5, specificity<0.01, min_frac=0.4, min_sites=30):
  A. final_score  (post-hoc from the base host_summary)
  B. specificity  (post-hoc from the base host_summary)
  C. min_sites / min_detect  (linkage-only re-runs via relink_sweep.py)
  D. min_frac (motif modified-fraction)  (linkage-only re-runs via relink_sweep.py)
Source Data written next to the figure.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE, RED, GREY = "#0072B2", "#D55E00", "#CC3311", "#888888"
S0, P0 = 0.5, 0.01
REPS = [f"orphan_300_rep{r}" for r in range(1, 6)]
sra = lambda x: str(x).split("_")[0]


def load_hs(lab, path):
    """load a host_summary.csv and tag orphan/correct using this rep's orphan list."""
    orph = set(l.strip() for l in open(f"{ROOT}/{lab}/{lab}.orphans.txt") if l.strip())
    h = pd.read_csv(path)
    h["is_orphan"] = h.MGE.isin(orph)
    h["correct"] = (~h.is_orphan) & (h.MGE.map(sra) == h.host.map(sra))
    return h, int((~h.is_orphan).sum()), int(h.is_orphan.sum())


def load_rep(lab):
    return load_hs(lab, f"{ROOT}/{lab}/modifi/{lab}/host_summary.csv")


def metrics(h, n_pl, n_or, s=S0, p=P0):
    """recall, PLANTED precision, orphan-FPR at operating point (s, p).

    precision counts only planted-ECE links (correct vs mislink); host-absent orphan
    false positives are the SEPARATE orphan-FPR axis, so the two error modes are not
    double-counted (matches the orphan-FPR figure's planted precision).
    """
    acc = h[(h.final_score > s) & (h.specificity < p)]
    tp = int(((~acc.is_orphan) & acc.correct).sum())
    planted_conf = int((~acc.is_orphan).sum())        # correct + mislink (planted only)
    orphFP = int(acc.is_orphan.sum())
    prec = tp / planted_conf if planted_conf else np.nan
    return tp / n_pl, prec, orphFP / n_or


def sweep_point(lab, tag):
    """(recall, precision, orphan-FPR) at the DEFAULT op point for one re-run host_summary."""
    hs = f"{ROOT}/{lab}/relink_sweep/{tag}/host_summary.csv"
    if not os.path.exists(hs):
        return (np.nan, np.nan, np.nan)
    return metrics(*load_hs(lab, hs))


def band(ax, x, mat, color, ls, label):
    """mat: reps x len(x); plot mean + 95% CI band."""
    m = np.nanmean(mat, 0)
    sem = np.nanstd(mat, 0, ddof=1) / np.sqrt(mat.shape[0])
    ax.fill_between(x, m - 1.96 * sem, m + 1.96 * sem, color=color, alpha=0.15, lw=0)
    ax.plot(x, m, color=color, lw=2, ls=ls, label=label)
    return m


def points(ax, x, mat, color, marker, ls, label):
    """mat: reps x len(x); plot mean markers + 95% CI error bars over discrete x."""
    m = np.nanmean(mat, 0)
    sem = np.nanstd(mat, 0, ddof=1) / np.sqrt(mat.shape[0])
    ax.errorbar(x, m, yerr=1.96 * sem, fmt=marker + ls, color=color, lw=2, ms=6,
                capsize=3, label=label)
    return m


MF_GRID = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]     # min_frac; default 0.4
MF_DEF = 0.4
MS_GRID = [10, 30, 60, 100, 200]             # min_sites (min_detect); default 30
MS_DEF = 30


def ms_tag(v):
    return "mf0.4" if v == MS_DEF else f"ms{v}"     # ms30 == default (mf0.4) re-run


def main():
    reps = [load_rep(l) for l in REPS]
    sgrid = np.linspace(0, 1, 81)
    pgrid = np.logspace(-4, 0, 81)

    # A/B: post-hoc from the base host_summary
    Rs = np.array([[metrics(h, npl, nor, s, P0) for s in sgrid] for h, npl, nor in reps])
    Rp = np.array([[metrics(h, npl, nor, S0, p) for p in pgrid] for h, npl, nor in reps])
    # C/D: from the linkage-only re-runs (score+spec at default in each)
    Rms = np.array([[sweep_point(lab, ms_tag(v)) for v in MS_GRID] for lab in REPS])
    Rmf = np.array([[sweep_point(lab, f"mf{v}") for v in MF_GRID] for lab in REPS])

    fig, axg = plt.subplots(2, 2, figsize=(13, 9.4))
    ax = axg.ravel()

    # A: final_score
    rec_s = band(ax[0], sgrid, Rs[:, :, 0], BLUE, "-", "recall")
    band(ax[0], sgrid, Rs[:, :, 1], ORANGE, "-", "precision")
    fpr_s = band(ax[0], sgrid, Rs[:, :, 2], RED, "--", "orphan-FPR")
    ax[0].axvline(S0, color=GREY, ls=":", lw=1)
    ax[0].text(S0, 1.03, "default 0.5", ha="center", va="bottom", fontsize=8, color=GREY)
    ax[0].set(xlabel="final_score accept threshold (> x)", ylabel="metric", ylim=(0, 1.12),
              title="a. Linkage score\n(specificity, min_frac, min_sites at default)")
    ax[0].legend(loc="center left", fontsize=9)

    # B: specificity
    band(ax[1], pgrid, Rp[:, :, 0], BLUE, "-", "recall")
    band(ax[1], pgrid, Rp[:, :, 1], ORANGE, "-", "precision")
    band(ax[1], pgrid, Rp[:, :, 2], RED, "--", "orphan-FPR")
    ax[1].axvline(P0, color=GREY, ls=":", lw=1)
    ax[1].text(P0, 1.03, "default 0.01", ha="center", va="bottom", fontsize=8, color=GREY)
    ax[1].set_xscale("log")
    ax[1].set(xlabel="specificity accept threshold (< x)", ylabel="metric", ylim=(0, 1.12),
              title="b. Specificity\n(final_score, min_frac, min_sites at default)")
    ax[1].legend(loc="center left", fontsize=9)

    # C: min_sites (min_detect), re-run
    points(ax[2], MS_GRID, Rms[:, :, 0], BLUE, "o", "-", "recall")
    points(ax[2], MS_GRID, Rms[:, :, 1], ORANGE, "s", "-", "precision")
    points(ax[2], MS_GRID, Rms[:, :, 2], RED, "^", "--", "orphan-FPR")
    ax[2].axvline(MS_DEF, color=GREY, ls=":", lw=1)
    ax[2].text(MS_DEF, 1.03, "default 30", ha="center", va="bottom", fontsize=8, color=GREY)
    ax[2].set(xlabel="min modified-site count per motif (min_sites)", ylabel="metric",
              ylim=(0, 1.12), title="c. Min modified-site count\n(final_score, specificity, min_frac at default)")
    ax[2].legend(loc="center left", fontsize=9)

    # D: min_frac (motif fraction), re-run
    points(ax[3], MF_GRID, Rmf[:, :, 0], BLUE, "o", "-", "recall")
    points(ax[3], MF_GRID, Rmf[:, :, 1], ORANGE, "s", "-", "precision")
    points(ax[3], MF_GRID, Rmf[:, :, 2], RED, "^", "--", "orphan-FPR")
    ax[3].axvline(MF_DEF, color=GREY, ls=":", lw=1)
    ax[3].text(MF_DEF, 1.03, "default 0.4", ha="center", va="bottom", fontsize=8, color=GREY)
    ax[3].set(xlabel="motif modified-fraction threshold (min_frac)", ylabel="metric",
              ylim=(0, 1.12), title="d. Motif fraction\n(final_score, specificity, min_sites at default)")
    ax[3].legend(loc="center left", fontsize=9)

    npl = np.mean([npl for _, npl, _ in reps]); nor = np.mean([nor for _, _, nor in reps])
    fig.suptitle(f"Threshold calibration with host-absent orphans as negatives "
                 f"(orphan_300, 5 reps; ~{npl:.0f} planted + ~{nor:.0f} orphans each; "
                 f"one panel per hyperparameter; mean +/- 95% CI)", fontsize=12.5, y=1.04)
    fig.tight_layout()
    out = f"{OUT}/orphan_calibration.pdf"

    # Source Data
    def cont_df(x, xname, M):
        n = M.shape[0]; d = {xname: x}
        for k, nm in enumerate(["recall", "precision", "orphan_fpr"]):
            d[f"{nm}_mean"] = np.nanmean(M[:, :, k], 0)
            d[f"{nm}_ci"] = 1.96 * np.nanstd(M[:, :, k], 0, ddof=1) / np.sqrt(n)
        return pd.DataFrame(d)
    cont_df(sgrid, "final_score", Rs).to_csv(f"{OUT}/orphan_calibration_score_sourcedata.csv", index=False)
    cont_df(pgrid, "specificity", Rp).to_csv(f"{OUT}/orphan_calibration_specificity_sourcedata.csv", index=False)
    cont_df(np.array(MS_GRID), "min_sites", Rms).to_csv(f"{OUT}/orphan_calibration_minsites_sourcedata.csv", index=False)
    cont_df(np.array(MF_GRID), "min_frac", Rmf).to_csv(f"{OUT}/orphan_calibration_minfrac_sourcedata.csv", index=False)

    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    di_s = np.argmin(np.abs(sgrid - S0))
    print(f"wrote {out}")
    print(f"at default: recall={rec_s[di_s]:.3f}  orphan-FPR={fpr_s[di_s]:.3f}")
    print(f"min_sites reps done: {np.sum(~np.isnan(Rms[:,:,0]),0)}")
    print(f"min_frac  reps done: {np.sum(~np.isnan(Rmf[:,:,0]),0)}")


if __name__ == "__main__":
    main()
