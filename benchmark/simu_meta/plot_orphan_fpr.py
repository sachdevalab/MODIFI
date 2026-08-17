#!/usr/bin/env python3
"""plot_orphan_fpr.py — C1 orphan-ECE false-positive analysis (host-removal negative control).

orphan_300 x5 reps: 29 donor species keep their host (planted ECEs) and 29 have the host
genome + chromosomal reads removed (orphan ECEs, true host absent), over 242 background.
Any confident link (final_score>0.5 & specificity<0.01) for an orphan is a false positive.
2x2 panels (mean +/- 95% CI over 5 reps unless noted):
  A. ROC: planted-correct-host score vs orphan best-host score (pooled), AUROC + default op point
  B. PR for the same detector, AUPRC
  C. per-rep metrics: orphan-FPR, planted recall, planted precision, overall FDR
  D. orphan-FPR by tier: novel-genus (no congener present) vs congener-present, Wilson CI
Source Data written next to the figure.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/orphan_fpr"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#BBBBBB"
S0, P0 = 0.5, 0.01
REPS = [f"orphan_300_rep{r}" for r in range(1, 6)]
sra = lambda x: str(x).split("_")[0]
genus = lambda s: (str(s).split()[0] if str(s) == str(s) and str(s) else "")


def mci(v):
    v = np.array([x for x in v if x == x], float)
    if not len(v):
        return np.nan, 0.0
    return v.mean(), (1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, 0.0, 0.0
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return p, max(0, c - half), min(1, c + half)


def load():
    rows, pos, neg, per_rep = [], [], [], []
    tier = {"novel_genus": [0, 0], "congener_present": [0, 0]}  # [fp, n]
    for lab in REPS:
        d = f"{ROOT}/{lab}"
        man = pd.read_csv(f"{d}/{lab}.manifest.csv")
        orph = set(x.strip() for x in open(f"{d}/{lab}.orphans.txt") if x.strip())
        h = pd.read_csv(f"{d}/modifi/{lab}/host_summary.csv")
        sp_of = man.drop_duplicates("sample").set_index("sample")["species"].to_dict()
        present = man[man.role.isin(["donor", "background"])]
        gn_present = set(genus(s) for s in present.species.dropna())

        h["pass"] = (h.final_score > S0) & (h.specificity < P0)
        h["orphan"] = h.MGE.isin(orph)
        h["correct"] = h.MGE.map(sra) == h.host.map(sra)
        pl = h[~h.orphan]
        orp = h[h.orphan].copy()
        tp = int((pl["pass"] & pl.correct).sum())
        ml = int((pl["pass"] & ~pl.correct).sum())
        ofp = int(orp["pass"].sum())
        rows.append(dict(rep=lab, planted_recall=tp / len(pl),
                         planted_prec=tp / (tp + ml) if (tp + ml) else np.nan,
                         orphan_fpr=ofp / len(orp),
                         fdr=(ml + ofp) / (tp + ml + ofp) if (tp + ml + ofp) else 0.0))
        rp = pl.loc[pl.correct, "final_score"].values
        rn = orp["final_score"].values
        pos += list(rp); neg += list(rn)
        per_rep.append((np.asarray(rp, float), np.asarray(rn, float)))
        orp["gp"] = orp.MGE.map(sra).map(sp_of).map(lambda s: genus(s) in gn_present)
        for _, rr in orp.iterrows():
            key = "congener_present" if rr.gp else "novel_genus"
            tier[key][1] += 1
            tier[key][0] += int(rr["pass"])
    return pd.DataFrame(rows), np.array(pos), np.array(neg), tier, per_rep


def roc_pr(pos, neg):
    thr = np.unique(np.concatenate([pos, neg, [1.01]]))[::-1]
    P, N = len(pos), len(neg)
    tpr, fpr, prec, rec = [], [], [], []
    for t in thr:
        tp = (pos >= t).sum(); fp = (neg >= t).sum()
        tpr.append(tp / P); fpr.append(fp / N)
        rec.append(tp / P); prec.append(tp / (tp + fp) if (tp + fp) else 1.0)
    fpr, tpr = np.array(fpr), np.array(tpr)
    rec, prec = np.array(rec), np.array(prec)
    trap = getattr(np, "trapezoid", getattr(np, "trapz", None))
    auroc = trap(tpr[np.argsort(fpr)], np.sort(fpr))
    o = np.argsort(rec)
    auprc = trap(prec[o], rec[o])
    return fpr, tpr, rec, prec, auroc, auprc, thr


def rep_roc(pos, neg, grid):
    """per-replicate ROC interpolated onto a common FPR grid; returns (tpr_on_grid, auroc)."""
    fpr, tpr, *_ , auroc, _, _ = roc_pr(pos, neg)
    o = np.argsort(fpr)
    it = np.interp(grid, fpr[o], tpr[o])
    it[0] = 0.0
    return it, auroc


def main():
    R, pos, neg, tier, per_rep = load()

    # per-replicate ROC on a common FPR grid -> mean curve + 95% CI band
    grid = np.linspace(0, 1, 201)
    tprs, aucs = [], []
    for rp, rn in per_rep:
        it, au = rep_roc(rp, rn, grid)
        tprs.append(it); aucs.append(au)
    tprs = np.vstack(tprs)
    mean_tpr = tprs.mean(0); mean_tpr[-1] = 1.0
    sem = tprs.std(0, ddof=1) / np.sqrt(tprs.shape[0])
    lo = np.clip(mean_tpr - 1.96 * sem, 0, 1)
    hi = np.clip(mean_tpr + 1.96 * sem, 0, 1)
    auroc_m, auroc_ci = mci(aucs)

    # default operating point per replicate (mean +/- 95% CI in both axes)
    op = np.array([[(rn > S0).mean(), (rp > S0).mean()] for rp, rn in per_rep])
    (ofpr_m, ofpr_ci), (orec_m, orec_ci) = mci(op[:, 0]), mci(op[:, 1])

    fig, axA = plt.subplots(figsize=(6, 5.4))
    for it in tprs:  # faint per-replicate curves
        axA.plot(grid, it, "-", color=BLUE, lw=0.8, alpha=0.25)
    axA.fill_between(grid, lo, hi, color=BLUE, alpha=0.20, lw=0,
                     label="95% CI across 5 reps")
    axA.plot(grid, mean_tpr, "-", color=BLUE, lw=2.2, label="mean ROC")
    axA.plot([0, 1], [0, 1], "--", color=GREY, lw=1)
    axA.set(xlabel="orphan false-positive rate", ylabel="planted true-host recall",
            xlim=(-0.02, 1.02), ylim=(0, 1.02),
            title=f"Real host link vs orphan (AUROC={auroc_m:.3f} +/- {auroc_ci:.3f})")
    axA.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    out = f"{OUT}/orphan_fpr_panel.pdf"

    # Source Data for the mean curve + band
    pd.DataFrame({"fpr_grid": grid, "mean_tpr": mean_tpr, "ci_lo": lo, "ci_hi": hi}) \
        .to_csv(f"{OUT}/orphan_roc_pr_sourcedata.csv", index=False)
    tk = ["novel_genus", "congener_present"]
    pd.DataFrame({"tier": tk, "orphan_fp": [tier[k][0] for k in tk],
                  "n_orphan": [tier[k][1] for k in tk],
                  "fpr": [tier[k][0] / tier[k][1] for k in tk]}) \
        .to_csv(f"{OUT}/orphan_tier_sourcedata.csv", index=False)

    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    print(f"AUROC={auroc_m:.3f}+/-{auroc_ci:.3f} | per-rep AUCs={[round(a,3) for a in aucs]}")
    print(f"default op: FPR={ofpr_m:.3f}+/-{ofpr_ci:.3f}, recall={orec_m:.3f}+/-{orec_ci:.3f}")
    for k in ("novel_genus", "congener_present"):
        print(f"  {k}: {tier[k][0]}/{tier[k][1]} = {tier[k][0]/tier[k][1]:.3f}")


if __name__ == "__main__":
    main()
