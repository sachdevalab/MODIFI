#!/usr/bin/env python3
"""plot_c3_coverage.py — ECE-host recall/precision vs donor sequencing depth (C3 titration).

Fixed bg_80 composition (40 donor species + 40 background) re-run at donor depth
{10,20,30,40x} (background fixed 10x), five replicates each (seeds 42-46). Locates the
usable-depth floor: recall rises with depth and saturates by ~30x while precision stays
~1.0. Mean +/- 95% CI across replicates; Source Data saved next to the figure.
"""
import os
import numpy as np
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


def rep_labels(base):
    labs = [base] + [f"{base}_rep{r}" for r in range(2, 6)]
    return [l for l in labs if os.path.exists(f"{ROOT}/{l}/modifi/{l}/host_summary.csv")]


def rep_metrics(lbl):
    h = pd.read_csv(f"{ROOT}/{lbl}/modifi/{lbl}/host_summary.csv")
    tot = pd.read_csv(f"{ROOT}/{lbl}/{lbl}.mge_list.tsv", sep="\t").seq_name.nunique()
    h["p"] = (h.final_score > 0.5) & (h.specificity < 0.01)
    tp = int((h["p"] & (h.MGE.map(sra) == h.host.map(sra))).sum())
    na = int(h["p"].sum())
    return tp / tot, (tp / na if na else np.nan)


def agg(vals):
    v = np.array([x for x in vals if x == x], float)
    if not len(v):
        return np.nan, 0.0
    return v.mean(), (1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)


def main():
    dep, rec_m, rec_ci, prec_m, prec_ci, nrep = [], [], [], [], [], []
    for d, base in POINTS:
        labs = rep_labels(base)
        if not labs:
            continue
        rr = [rep_metrics(l) for l in labs]
        rm, rc = agg([x[0] for x in rr]); pm, pc = agg([x[1] for x in rr])
        dep.append(d); rec_m.append(rm); rec_ci.append(rc)
        prec_m.append(pm); prec_ci.append(pc); nrep.append(len(labs))

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.errorbar(dep, rec_m, yerr=rec_ci, fmt="-o", color=BLUE, lw=2, ms=8, capsize=4, label="recall")
    ax.errorbar(dep, prec_m, yerr=prec_ci, fmt="-s", color=ORANGE, lw=2, ms=7, capsize=4, label="precision")
    for x, r, p in zip(dep, rec_m, prec_m):
        ax.annotate(f"{r:.2f}", (x, r), textcoords="offset points", xytext=(0, -16), ha="center", fontsize=9, color=BLUE)
        ax.annotate(f"{p:.2f}", (x, p), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=9, color=ORANGE)
    ax.set(xlabel="donor sequencing depth (x)", ylabel="metric", ylim=(0, 1.05),
           title=f"C3: ECE-host recall & precision vs coverage\n(bg_80 composition; background fixed 10x; "
                 f"mean +/- 95% CI, n={min(nrep)}-{max(nrep)} reps)")
    ax.set_xticks(dep)
    ax.legend(loc="center right", fontsize=10)
    fig.tight_layout()
    out = f"{OUT}/c3_coverage_titration.pdf"

    pd.DataFrame({"donor_depth_x": dep, "n_rep": nrep,
                  "recall_mean": rec_m, "recall_ci": rec_ci,
                  "precision_mean": prec_m, "precision_ci": prec_ci}) \
        .to_csv(out.replace(".pdf", "_sourcedata.csv"), index=False)
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out} | reps: {dict(zip(dep, nrep))}")
    print(f"depth->recall: {dict(zip(dep, [round(r,2) for r in rec_m]))}")


if __name__ == "__main__":
    main()
