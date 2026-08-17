#!/usr/bin/env python3
"""plot_benchmark_robustness.py — integrated comprehensive benchmark of MODIFI ECE-host
linkage across three axes, built by mixing real isolate CCS reads (real IPD kinetics, exact
SRA-prefix ground truth; no simulated signal):
  a. community complexity   (C1): 25 -> 300 genomes, reference mode, species-level
  b. sequencing depth       (C3): donor 10 -> 40x, fixed bg_80 composition
  c. end-to-end de-novo      (C4): mix -> hifiasm_meta -> ECE calling -> MODIFI, recall over ALL ECEs
Each point is mean +/- 95% CI across 5 replicates (seeds 42-46). Source Data saved alongside.
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C1 = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
E2E = "/home/shuaiw/borg/paper/simu_meta_dir"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/benchmark_summary"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE = "#0072B2", "#D55E00"
sra = lambda x: str(x).split("_")[0]


def agg(v):
    v = np.array([x for x in v if x == x], float)
    if not len(v):
        return np.nan, 0.0
    return v.mean(), (1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)


def reps_c1(base):
    L = [base] + [f"{base}_rep{r}" for r in range(2, 6)]
    return [l for l in L if os.path.exists(f"{C1}/{l}/modifi/{l}/host_summary.csv")]


def metric_c1(lbl):
    h = pd.read_csv(f"{C1}/{lbl}/modifi/{lbl}/host_summary.csv")
    tot = pd.read_csv(f"{C1}/{lbl}/{lbl}.mge_list.tsv", sep="\t").seq_name.nunique()
    p = (h.final_score > 0.5) & (h.specificity < 0.01)
    tp = int((p & (h.MGE.map(sra) == h.host.map(sra))).sum()); na = int(p.sum())
    return tp / tot, (tp / na if na else np.nan)


def series_c1(points):
    xs, rm, rc, pm, pc, nr = [], [], [], [], [], []
    for base, x in points:
        R = reps_c1(base)
        if not R:
            continue
        rr = [metric_c1(l) for l in R]
        a, b = agg([v[0] for v in rr]); c, d = agg([v[1] for v in rr])
        xs.append(x); rm.append(a); rc.append(b); pm.append(c); pc.append(d); nr.append(len(R))
    return xs, rm, rc, pm, pc, nr


def _parse_ece(cell):
    return [x for x in str(cell).replace(";", ",").split(",") if x.strip()] if isinstance(cell, str) and cell.strip() else []


def series_e2e(points):
    xs, rm, rc, pm, pc, nr = [], [], [], [], [], []
    for base, x in points:
        dirs = [d for d in sorted(glob.glob(f"{E2E}/C4_{base}") + glob.glob(f"{E2E}/C4_{base}_rep*"))
                if os.path.exists(f"{d}/e2e_predictions.csv")]
        if not dirs:
            continue
        rr, pp = [], []
        for d in dirs:
            man = pd.read_csv(f"{d}/toy.manifest.csv")
            ecol = [c for c in man.columns if "ece" in c.lower() and "contig" in c.lower()]
            all_ece = set()
            for v in man[ecol[0]]:
                all_ece.update(_parse_ece(v))
            et = pd.read_csv(f"{d}/ece_truth.tsv", sep="\t")
            a2c = et.set_index("ece_contig")["matched_curated_ece"].to_dict()
            x2 = pd.read_csv(f"{d}/e2e_predictions.csv"); x2["cur"] = x2.MGE.map(a2c)
            rr.append(x2[x2.correct].cur.dropna().nunique() / len(all_ece) if all_ece else np.nan)
            pp.append(x2.correct.sum() / x2.assigned.sum() if x2.assigned.sum() else np.nan)
        a, b = agg(rr); c, d = agg(pp)
        xs.append(x); rm.append(a); rc.append(b); pm.append(c); pc.append(d); nr.append(len(dirs))
    return xs, rm, rc, pm, pc, nr


def panel(ax, xs, rm, rc, pm, pc, xticklabels, xlabel, title, logx=False):
    ax.errorbar(xs, rm, yerr=rc, fmt="-o", color=BLUE, lw=2, ms=7, capsize=4, label="recall")
    ax.errorbar(xs, pm, yerr=pc, fmt="-s", color=ORANGE, lw=2, ms=6, capsize=4, label="precision")
    if logx:
        ax.set_xscale("log")
    ax.set_xticks(xs); ax.set_xticklabels(xticklabels)
    ax.set(ylim=(0, 1.08), xlabel=xlabel, ylabel="metric", title=title)
    ax.legend(loc="lower left", fontsize=9)


def main():
    C1_PTS = [("ladder_25", 25), ("ladder_40", 40), ("ladder_58", 58),
              ("bg_80", 80), ("bg_150", 150), ("bg_300", 300)]
    COV_PTS = [("cov_d10", 10), ("cov_d20", 20), ("bg_80", 30), ("cov_d40", 40)]
    E2E_PTS = [("bg_80", 80), ("bg_150", 150), ("bg_300", 300)]

    ca = series_c1(C1_PTS)
    cb = series_c1(COV_PTS)
    cc = series_e2e(E2E_PTS)

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))
    panel(ax[0], ca[0], ca[1], ca[2], ca[3], ca[4], [str(x) for x in ca[0]],
          "community size (genomes)", "a. Community complexity\n(reference mode, species level)", logx=True)
    ax[0].axvspan(24, 58, color="#EEEEEE", alpha=0.6, zorder=0)
    panel(ax[1], cb[0], cb[1], cb[2], cb[3], cb[4], [str(x) for x in cb[0]],
          "donor sequencing depth (x)", "b. Sequencing depth\n(fixed bg_80 composition, background 10x)")
    panel(ax[2], cc[0], cc[1], cc[2], cc[3], cc[4], [str(x) for x in cc[0]],
          "community size (genomes)", "c. End-to-end de-novo\n(assembly + ECE calling; recall over all ECEs)")

    nmin = min([min(s[5]) for s in (ca, cb, cc) if s[5]]); nmax = max([max(s[5]) for s in (ca, cb, cc) if s[5]])
    fig.suptitle(f"Comprehensive benchmark of MODIFI ECE-host linkage (real isolate reads, exact ground truth; "
                 f"mean +/- 95% CI, n={nmin}-{nmax} replicates)", fontsize=12.5, y=1.03)
    fig.tight_layout()
    out = f"{OUT}/benchmark_robustness.pdf"

    rows = []
    for name, s, xl in [("complexity", ca, "n_genomes"), ("coverage", cb, "donor_depth_x"), ("end2end_denovo", cc, "n_genomes")]:
        for i, x in enumerate(s[0]):
            rows.append(dict(axis=name, x=x, x_label=xl, n_rep=s[5][i],
                             recall_mean=s[1][i], recall_ci=s[2][i],
                             precision_mean=s[3][i], precision_ci=s[4][i]))
    pd.DataFrame(rows).to_csv(out.replace(".pdf", "_sourcedata.csv"), index=False)
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    for r in rows:
        print(f"  {r['axis']:16s} x={r['x']:>4} n={r['n_rep']} recall={r['recall_mean']:.2f} prec={r['precision_mean']:.2f}")


if __name__ == "__main__":
    main()
