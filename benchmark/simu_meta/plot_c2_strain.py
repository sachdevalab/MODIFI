#!/usr/bin/env python3
"""plot_c2_strain.py — C2 strain-mixture: recall & precision vs con-specific strain depth K.

Fixed 58-donor-species + 242-background community; K = max con-specific strains/species
{1,2,3,4,all}. Each K adds more con-specific strains to the same 58 donor species; the
x-axis is annotated with the resulting number of donor strains. Metrics at the default
operating point (score>0.5 & spec<0.01). rep1.
"""
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SS = PRED = None
for root, _, files in os.walk("/home/shuaiw/borg/paper/simu_meta_dir/C1"):
    if "/eval" in root:
        if "size_summary.csv" in files:
            SS = os.path.join(root, "size_summary.csv")
        if "predictions.csv" in files:
            PRED = os.path.join(root, "predictions.csv")
ROOT = "/home/shuaiw/borg/paper/simu_meta_dir/C1"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/strain_het"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
KMAP = [("bg_300", "1", 1), ("bg300_k2", "2", 2), ("bg300_k3", "3", 3),
        ("bg300_k4", "4", 4), ("bg300_kall", "all", 5)]
BLUE, ORANGE = "#0072B2", "#D55E00"


def strain_counts(lbl):
    """(donor strains, total genomes) in the community."""
    m = f"{ROOT}/{lbl}/{lbl}.manifest.csv"
    if not os.path.exists(m):
        return 0, 0
    d = pd.read_csv(m)
    return int((d["role"] == "donor").sum()), int(len(d))


def n_donor_strains(lbl):
    return strain_counts(lbl)[0]


def metrics(pred, labels):
    """Per K, at the operating point: species-level and strain-level recall/precision, and
    strain accuracy (among confident correct-species calls, fraction to the exact strain)."""
    p = pred[pred.setting.isin(labels)].copy()
    conf = p[(p.final_score > 0.5) & (p.specificity < 0.01)]
    out = {}
    for lbl in labels:
        n = int((p.setting == lbl).sum())                      # planted ECEs (recall denom)
        g = conf[conf.setting == lbl]
        ci = int((g.resolution == "correct_isolate").sum())     # exact strain
        ws = int((g.resolution == "correct_species_wrong_strain").sum())
        assigned = int(len(g))
        sp_tp = ci + ws                                         # species-correct (either strain)
        out[lbl] = dict(
            sp_recall=sp_tp / n, sp_prec=sp_tp / assigned if assigned else float("nan"),
            st_recall=ci / n,    st_prec=ci / assigned if assigned else float("nan"),
            strain_acc=100 * ci / (ci + ws) if (ci + ws) else float("nan"))
    return out


def _rp(ax, x, rec, prec, title, ylab):
    ax.plot(x, rec, "-o", color=BLUE, lw=2, ms=8, label="recall")
    ax.plot(x, prec, "-s", color=ORANGE, lw=2, ms=7, label="precision")
    for xi, r, p in zip(x, rec, prec):
        ax.annotate(f"{r:.2f}", (xi, r), textcoords="offset points", xytext=(0, -16), ha="center", fontsize=8.5, color=BLUE)
        ax.annotate(f"{p:.2f}", (xi, p), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8.5, color=ORANGE)
    ax.set(ylim=(0, 1.05), title=title, ylabel=ylab)
    ax.legend(loc="lower left", fontsize=9)


def main():
    pred = pd.read_csv(PRED)
    K, x, nstr, labels = [], [], [], []
    for lbl, k, o in KMAP:
        m = f"{ROOT}/{lbl}/{lbl}.manifest.csv"
        if not os.path.exists(m):
            continue
        K.append(k); x.append(o); nstr.append(n_donor_strains(lbl)); labels.append(lbl)
    M = metrics(pred, labels)
    ndon = nstr
    ntot = [strain_counts(l)[1] for l in labels]
    xt = [f"K={k}" for k in K]
    xl = "max con-specific strains per donor species (K)"

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    axA, axB, axC, axD = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]

    # A. strain-count bar plot (donor strains + total genomes) per K
    import numpy as np
    xi = np.arange(len(x))
    axA.bar(xi - 0.2, ndon, 0.38, color="#0072B2", label="donor strains")
    axA.bar(xi + 0.2, ntot, 0.38, color="#BBBBBB", label="total genomes (incl. 242 background)")
    for j, (d, t) in enumerate(zip(ndon, ntot)):
        axA.annotate(str(d), (j - 0.2, d), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8)
        axA.annotate(str(t), (j + 0.2, t), textcoords="offset points", xytext=(0, 3), ha="center", fontsize=8)
    axA.set_xticks(xi); axA.set_xticklabels(xt)
    axA.set(xlabel=xl, ylabel="count", title="A. Community strain composition vs K")
    axA.legend(loc="upper left", fontsize=8.5)

    _rp(axB, x, [M[l]["sp_recall"] for l in labels], [M[l]["sp_prec"] for l in labels],
        "B. Species-level recall & precision", "species-level metric")
    _rp(axC, x, [M[l]["st_recall"] for l in labels], [M[l]["st_prec"] for l in labels],
        "C. Strain-level recall & precision", "strain-level metric")
    sa = [M[l]["strain_acc"] for l in labels]
    axD.plot(x, sa, "-o", color="#009E73", lw=2, ms=8)
    for xv, y in zip(x, sa):
        axD.annotate(f"{y:.1f}%", (xv, y), textcoords="offset points", xytext=(0, 9), ha="center", fontsize=9)
    axD.set(ylim=(80, 101), title="D. Strain accuracy (of correct-species calls)", ylabel="strain accuracy (%)")
    for a in (axB, axC, axD):
        a.set_xticks(x); a.set_xticklabels(xt); a.set_xlabel(xl)
    fig.suptitle("C2: strain-mixture ECE-host assignment (58 donor species + 242 background; rep1)",
                 fontsize=12.5, y=1.02)
    fig.tight_layout()
    out = f"{OUT}/c2_strain_recall_precision.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    # remove the old (resolution-stacked) figure name so the folder isn't ambiguous
    for e in (".pdf", ".png"):
        old = f"{OUT}/c2_strain_resolution{e}"
        if os.path.exists(old):
            os.remove(old)
    print(f"wrote {out}  | donor strains per K: {dict(zip(K, nstr))}")


if __name__ == "__main__":
    main()
