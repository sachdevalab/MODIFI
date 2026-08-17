#!/usr/bin/env python3
"""plot_e2e_robustness.py — MODIFI is robust end-to-end (de-novo pipeline).

Runs the realistic path (mix real reads -> hifiasm_meta co-assembly -> ground-truth ECE
selection -> MODIFI contig-level linkage) on the bg_80/150/300 communities, pooling all
available replicates. Shows that even through de-novo assembly and ECE calling MODIFI keeps
ECE->host precision near 1.0 with stable recall as community complexity grows from 80 to 300
genomes. Recall is scored over ALL curated ECEs present in the community (assembly-lost ECEs
count as misses, the honest end-to-end denominator); precision is over the confident assembled
links. An ECE is correct when its best host contig traces to the ECE's own source isolate at
the default operating point (final_score > 0.5, specificity < 0.01). Source Data saved alongside.
"""
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/end2end"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "#BBBBBB"
COMMS = [("bg_80", 80), ("bg_150", 150), ("bg_300", 300)]     # label, total genomes


def rep_dirs(base):
    ds = sorted(glob.glob(f"{ROOT}/C4_{base}") + glob.glob(f"{ROOT}/C4_{base}_rep*"))
    return [d for d in ds if os.path.exists(f"{d}/e2e_predictions.csv")]


def _parse_ece(cell):
    if not isinstance(cell, str) or not cell.strip():
        return []
    return [x for x in cell.replace(";", ",").split(",") if x.strip()]


def rep_metrics(d):
    """Recall is over ALL curated ECEs present in the community (assembly-lost ECEs count
    as misses); precision is over the confident assembled links (unaffected by assembly loss)."""
    x = pd.read_csv(f"{d}/e2e_predictions.csv")
    man = pd.read_csv(f"{d}/toy.manifest.csv")
    ecol = [c for c in man.columns if "ece" in c.lower() and "contig" in c.lower()]
    all_ece = set()
    if ecol:
        for v in man[ecol[0]]:
            all_ece.update(_parse_ece(v))
    et = pd.read_csv(f"{d}/ece_truth.tsv", sep="\t")
    a2c = et.set_index("ece_contig")["matched_curated_ece"].to_dict()   # assembled contig -> curated ECE
    x["curated"] = x.MGE.map(a2c)

    n_all = len(all_ece)
    correct_curated = int(x[x.correct].curated.dropna().nunique())      # distinct curated ECEs correctly linked
    conf = int(x.assigned.sum()); corr_contig = int(x.correct.sum())
    return dict(n_all=n_all, n_assembled_ece=x.curated.nunique(),
                recall=correct_curated / n_all if n_all else np.nan,     # over ALL ECEs
                precision=corr_contig / conf if conf else np.nan,        # over confident assembled links
                correct_curated=correct_curated, asm_recovery=x.curated.nunique() / n_all if n_all else np.nan)


def agg(v):
    v = np.array([x for x in v if x == x], float)
    if not len(v):
        return np.nan, 0.0
    return v.mean(), (1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0)


def main():
    sizes, xs, nrep = [], [], []
    rec_m, rec_ci, prec_m, prec_ci = [], [], [], []
    src_rows = []
    for i, (base, size) in enumerate(COMMS):
        dirs = rep_dirs(base)
        if not dirs:
            continue
        ms = [rep_metrics(d) for d in dirs]
        rm, rc = agg([m["recall"] for m in ms]); pm, pc = agg([m["precision"] for m in ms])
        am, _ = agg([m["asm_recovery"] for m in ms])
        sizes.append(size); xs.append(i); nrep.append(len(dirs))
        rec_m.append(rm); rec_ci.append(rc); prec_m.append(pm); prec_ci.append(pc)
        src_rows.append(dict(community=base, total_genomes=size, n_rep=len(dirs),
                             n_all_curated_ECE=int(np.mean([m["n_all"] for m in ms])),
                             recall_all_ECE_mean=rm, recall_all_ECE_ci=rc,
                             precision_mean=pm, precision_ci=pc, assembly_recovery_mean=am))

    fig, ax0 = plt.subplots(figsize=(6.6, 5))

    # recall (over ALL curated ECEs) & precision vs community size
    ax0.errorbar(xs, rec_m, yerr=rec_ci, fmt="-o", color=BLUE, lw=2, ms=8, capsize=4, label="recall (all ECEs)")
    ax0.errorbar(xs, prec_m, yerr=prec_ci, fmt="-s", color=ORANGE, lw=2, ms=7, capsize=4, label="precision")
    for x, r, p in zip(xs, rec_m, prec_m):
        ax0.annotate(f"{r:.2f}", (x, r), textcoords="offset points", xytext=(0, -16), ha="center", fontsize=9, color=BLUE)
        ax0.annotate(f"{p:.2f}", (x, p), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=9, color=ORANGE)
    ax0.set_xticks(xs); ax0.set_xticklabels([f"{s}\ngenomes" for s in sizes])
    ax0.set(ylim=(0, 1.08), xlabel="community size", ylabel="metric",
            title=f"MODIFI de-novo ECE-host linkage is robust to complexity\n"
                  f"(end-to-end pipeline; recall over all ECEs; mean +/- 95% CI, n={min(nrep)}-{max(nrep)} reps)")
    ax0.legend(loc="center left", fontsize=10)
    fig.tight_layout()
    out = f"{OUT}/e2e_robustness.pdf"
    pd.DataFrame(src_rows).to_csv(out.replace(".pdf", "_sourcedata.csv"), index=False)
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out}")
    for r in src_rows:
        print(f"  {r['community']:6s} n_rep={r['n_rep']} nECE={r['n_all_curated_ECE']} "
              f"recall_all={r['recall_all_ECE_mean']:.2f} prec={r['precision_mean']:.2f} "
              f"asm_recovery={r['assembly_recovery_mean']:.2f}")


if __name__ == "__main__":
    main()
