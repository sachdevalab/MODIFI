#!/usr/bin/env python3
"""plot_c2_strain.py — C2 strain-mixture: ECE-host assignment vs con-specific strain depth K.

Fixed 58-donor-species + 242-background community; K = max con-specific strains/species
{1,2,3,4,all} (K=1 = bg_300). Metrics at the default operating point (score>0.5 & spec<0.01),
per replicate, shown as mean +/- 95% CI. 2x2 panels:
  A. community strain composition (donor strains + total genomes) per K
  B. species-level recall & precision
  C. strain-level recall & precision
  D. strain accuracy (of confident correct-species calls, fraction to the exact strain)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import build_community as bc

ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/strain_het"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE = "#0072B2", "#D55E00"
S0, P0 = 0.5, 0.01
# K -> ordered rep community labels (rep1 = bare; rep2-5 = _repN)
KDEF = [("1", "strain_mix_k1"), ("2", "strain_mix_k2"), ("3", "strain_mix_k3"),
        ("4", "strain_mix_k4"), ("all", "strain_mix_kall")]


def sra(x):
    return str(x).split("_")[0]


def rep_labels(base):
    labs = [f"{base}_rep{r}" for r in range(1, 6)]
    return [l for l in labs if os.path.exists(f"{ROOT}/{l}/modifi/{l}/host_summary.csv")]


def rep_metrics(label, sp_of):
    """operating-point metrics for one replicate community."""
    d = f"{ROOT}/{label}"
    h = pd.read_csv(f"{d}/modifi/{label}/host_summary.csv")
    tot = pd.read_csv(f"{d}/{label}.mge_list.tsv", sep="\t")["seq_name"].nunique()
    h["pass"] = (h.final_score > S0) & (h.specificity < P0)
    ei, hi = h.MGE.map(sra), h.host.map(sra)
    cs = h["pass"] & (ei == hi)                                   # correct strain (exact isolate)
    csp = h["pass"] & (ei.map(sp_of) == hi.map(sp_of)) & ei.map(sp_of).notna()  # correct species
    na = int(h["pass"].sum())
    return dict(sp_rec=csp.sum() / tot, sp_prec=csp.sum() / na if na else np.nan,
                st_rec=cs.sum() / tot, st_prec=cs.sum() / na if na else np.nan,
                strain_acc=100 * cs.sum() / csp.sum() if csp.sum() else np.nan)


def agg(vals):
    v = np.array([x for x in vals if x == x], float)
    if len(v) == 0:
        return np.nan, 0.0
    ci = 1.96 * v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0.0
    return v.mean(), ci


def strain_counts(base):
    m = pd.read_csv(f"{ROOT}/{base}_rep1/{base}_rep1.manifest.csv")
    return int((m.role == "donor").sum()), int(len(m))


def main():
    donors, bg, _ = bc.load_pool()
    sp_of = pd.concat([donors[["Sample", "Species"]], bg[["Sample", "Species"]]]) \
        .drop_duplicates("Sample").set_index("Sample")["Species"].to_dict()

    K, x, ndon, ntot, nrep = [], [], [], [], []
    M = {k: {"sp_rec": [], "sp_prec": [], "st_rec": [], "st_prec": [], "strain_acc": []} for k, _ in KDEF}
    for i, (k, base) in enumerate(KDEF):
        labs = rep_labels(base)
        if not labs:
            continue
        for lab in labs:
            m = rep_metrics(lab, sp_of)
            for key in M[k]:
                M[k][key].append(m[key])
        K.append(k); x.append(i); nrep.append(len(labs))
        dn, tn = strain_counts(base); ndon.append(dn); ntot.append(tn)

    def band(ax, key, color, marker, label):
        means, cis = zip(*[agg(M[k][key]) for k in K])
        means, cis = np.array(means), np.array(cis)
        ax.errorbar(x, means, yerr=cis, fmt=marker + "-", color=color, lw=2, ms=7,
                    capsize=3, label=label)
        return means

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    axA, axB, axC, axD = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]
    xt = [f"K={k}" for k in K]; xl = "max con-specific strains per donor species (K)"

    # A: strain composition
    xi = np.arange(len(x))
    axA.bar(xi - 0.2, ndon, 0.38, color=BLUE, label="donor strains")
    axA.bar(xi + 0.2, ntot, 0.38, color="#BBBBBB", label="total genomes (incl. 242 bg)")
    for j, (dd, tt) in enumerate(zip(ndon, ntot)):
        axA.annotate(str(dd), (j - 0.2, dd), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8)
        axA.annotate(str(tt), (j + 0.2, tt), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8)
    axA.set_xticks(xi); axA.set_xticklabels(xt)
    axA.set(xlabel=xl, ylabel="count", title="A. Community strain composition vs K")
    axA.legend(loc="upper left", fontsize=8.5)

    band(axB, "sp_rec", BLUE, "o", "recall"); band(axB, "sp_prec", ORANGE, "s", "precision")
    axB.set(ylim=(0, 1.05), title="B. Species-level recall & precision", ylabel="species-level metric")
    axB.legend(loc="lower left", fontsize=9)
    band(axC, "st_rec", BLUE, "o", "recall"); band(axC, "st_prec", ORANGE, "s", "precision")
    axC.set(ylim=(0, 1.05), title="C. Strain-level recall & precision", ylabel="strain-level metric")
    axC.legend(loc="lower left", fontsize=9)
    sa = band(axD, "strain_acc", "#009E73", "o", "strain accuracy")
    for xv, y in zip(x, sa):
        if y == y:
            axD.annotate(f"{y:.1f}%", (xv, y), textcoords="offset points", xytext=(0, -15),
                         ha="center", va="top", fontsize=9)
    axD.set(ylim=(0, 100), title="D. Strain accuracy (of correct-species calls)", ylabel="strain accuracy (%)")
    for a in (axB, axC, axD):
        a.set_xticks(x); a.set_xticklabels(xt); a.set_xlabel(xl)

    fig.suptitle(f"C2: strain-mixture ECE-host assignment (58 donor species + 242 background, clean split; "
                 f"mean +/- 95% CI, n={min(nrep)}-{max(nrep)} reps)", fontsize=12.5, y=1.02)
    fig.tight_layout()
    out = f"{OUT}/strain_mix_recall_precision.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out} | reps per K: {dict(zip(K, nrep))}")


if __name__ == "__main__":
    main()
