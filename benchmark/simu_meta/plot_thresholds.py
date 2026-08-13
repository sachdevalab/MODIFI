#!/usr/bin/env python3
"""
plot_thresholds.py — ROC calibration of MODIFI's linkage thresholds (Part G / C5)
across ALL completed ladder + background communities.

One ROC line per community setting (ladder_10 ... bg_300), coloured by community size;
the 5 replicates (seeds 42-46) of each setting are averaged onto a common FPR grid and
shown as mean +/- 1 SD shaded band. Ground truth: an ECE-host pair is correct iff same
SRA prefix. Two panels: `final_score` and `specificity` as the discriminating parameter.

Vector PDF -> /home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold/
Re-run any time; bands fill in as more replicates finish.
"""
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from sklearn.metrics import roc_curve, auc

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3})

# settings to include, in size order (label -> approx genome count for colouring)
SIZE = {"ladder_10": 10, "ladder_25": 25, "ladder_40": 40, "ladder_58": 58,
        "bg_80": 80, "bg_150": 153, "bg_300": 304}
MEAN_FPR = np.linspace(0, 1, 300)


def sra(x): return str(x).split("_")[0]


def find_reps(setting):
    """All completed replicate dirs for a setting: rep1 = bare label, reps>=2 = _repN."""
    reps = []
    for d in sorted(glob.glob(f"{ROOT}/{setting}") + glob.glob(f"{ROOT}/{setting}_rep*")):
        label = os.path.basename(d)
        if label != setting and not re.match(rf"^{setting}_rep\d+$", label):
            continue
        hs = f"{d}/modifi/{label}/host_summary.csv"
        hosts = f"{d}/modifi/{label}/hosts"
        if os.path.isdir(hosts):
            reps.append((label, hosts))
    return reps


def mean_ece_frac(motif_info):
    """Mean ECE (plasmid) methylation fraction across a pair's shared motifs, from
    motif_info (motif:cp:host_total:host_meth:pl_total:pl_meth:occ:occ_len)."""
    if not isinstance(motif_info, str) or not motif_info:
        return 0.0
    fr = []
    for e in motif_info.split(";"):
        p = e.split(":")
        if len(p) >= 6:
            try:
                pt, pm = int(p[4]), int(p[5])
                if pt > 0:
                    fr.append(pm / pt)
            except ValueError:
                pass
    return float(np.mean(fr)) if fr else 0.0


def rep_pairs(hosts_dir):
    rows = []
    for f in glob.glob(f"{hosts_dir}/*.host_prediction.csv"):
        try:
            d = pd.read_csv(f)
        except Exception:
            continue
        if d.empty:
            continue
        d["correct"] = d["MGE"].map(sra) == d["host"].map(sra)
        d["ece_frac"] = d["motif_info"].apply(mean_ece_frac)
        rows.append(d[["MGE", "host", "final_score", "specificity",
                       "total_sites", "ece_frac", "correct"]])
    return pd.concat(rows, ignore_index=True) if rows else None


# param -> (column, higher_is_positive transform)
def _score_vec(df, param):
    if param == "final_score":
        return df["final_score"].values.astype(float)
    if param == "specificity":  # lower specificity => more likely true
        return -np.log10(np.clip(df["specificity"].values.astype(float), 1e-6, None))
    if param == "total_sites":
        return df["total_sites"].values.astype(float)
    if param == "ece_frac":
        return df["ece_frac"].values.astype(float)
    raise ValueError(param)


def roc_for(df, param):
    y = df["correct"].values.astype(bool)
    if y.sum() == 0 or (~y).sum() == 0:
        return None
    fpr, tpr, _ = roc_curve(y, _score_vec(df, param))
    interp = np.interp(MEAN_FPR, fpr, tpr); interp[0] = 0.0
    return interp, auc(fpr, tpr)


def panel(ax, param, title):
    norm = colors.LogNorm(vmin=10, vmax=304)
    cmap = cm.viridis
    for setting in [s for s in SIZE if find_reps(s)]:
        tprs, aucs = [], []
        for label, hosts in find_reps(setting):
            df = rep_pairs(hosts)
            if df is None:
                continue
            r = roc_for(df, param)
            if r is None:
                continue
            tprs.append(r[0]); aucs.append(r[1])
        if not tprs:
            continue
        tprs = np.vstack(tprs)
        mean_tpr = tprs.mean(0)
        col = cmap(norm(SIZE[setting]))
        nrep = len(tprs)
        lab = f"{setting} ({SIZE[setting]}g) AUC={np.mean(aucs):.3f}"
        if nrep > 1:
            lab += f"±{np.std(aucs):.3f} (n={nrep})"
        ax.plot(MEAN_FPR, mean_tpr, color=col, lw=1.8, label=lab)
        if nrep > 1:
            sd = tprs.std(0)
            ax.fill_between(MEAN_FPR, np.clip(mean_tpr - sd, 0, 1),
                            np.clip(mean_tpr + sd, 0, 1), color=col, alpha=0.18, lw=0)
    ax.plot([0, 1], [0, 1], "--", color="#999", lw=1)
    ax.set(xlabel="False-positive rate", ylabel="True-positive rate (recall)",
           title=title, xlim=(-0.01, 1.0), ylim=(0, 1.01))
    ax.legend(loc="lower right", fontsize=7.5)


def main():
    have = [s for s in SIZE if find_reps(s)]
    print(f"[plot] settings with data: {have}")
    panels = [
        ("final_score", "A. linkage score  (cited: >0.8)"),
        ("specificity", "B. specificity  (cited: <0.001)"),
        ("ece_frac", "C. motif methylation fraction  (cited: >0.3)"),
        ("total_sites", "D. # modified sites  (cited: ≥100)"),
    ]
    fig, ax = plt.subplots(2, 2, figsize=(13, 11))
    for a, (param, title) in zip(ax.ravel(), panels):
        panel(a, param, "ROC — " + title)
        a.set_xlim(-0.005, 0.30)   # zoom: curves hug the corner
    fig.suptitle("MODIFI linkage-threshold discrimination across community complexity "
                 "(all 4 reviewer parameters; per community, band = ±SD over replicates)",
                 fontsize=12.5, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    pdf = f"{OUT}/threshold_roc.pdf"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(pdf.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"[plot] wrote {pdf}")


if __name__ == "__main__":
    main()
