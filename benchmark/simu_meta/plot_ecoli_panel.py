#!/usr/bin/env python3
"""plot_ecoli_panel.py — strain-mixture assembly stress test (E. coli panel).

Panel A: de-novo assembly of 8 con-specific E. coli strains produces inter-strain chimeras
(most E. coli contigs), yet CheckM2 still reports ~100% completeness (masking the problem).
Panel B: despite the chimeric assembly, MODIFI keeps precision 1.0 for strain-level ECE
assignment (reference vs de-novo); recall is low because con-specific strains share most
motifs, so ambiguous ECEs are declined, not mis-assigned.
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REF = "/home/shuaiw/borg/paper/simu_meta_dir/C1/ecoli_panel"
E2E = "/home/shuaiw/borg/paper/simu_meta_dir/C4_ecoli_panel"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/e_coli"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
BLUE, ORANGE, GREY, RED = "#0072B2", "#E69F00", "#BBBBBB", "#D55E00"


def sra(x):
    return str(x).split("_")[0]


def main():
    man = pd.read_csv(f"{REF}/ecoli_panel.manifest.csv")
    ecoli = set(man[man.role == "donor"]["sample"])
    # assembly: fraction of E. coli contigs that are chimeric; CheckM2 completeness
    org = pd.read_csv(f"{E2E}/contig_origin.tsv", sep="\t")
    e = org[org.origin_sample.isin(ecoli)]
    chim_pct = 100 * e.chimeric.mean()
    cm = pd.read_csv(f"{E2E}/checkm2/quality_report.tsv", sep="\t")
    comp = cm[cm.Name.isin(ecoli)].Completeness.median()

    # de-novo strain-level linkage
    d = pd.read_csv(f"{E2E}/e2e_predictions.csv")
    d["p"] = (d.final_score > 0.5) & (d.specificity < 0.01)
    ok = d["p"] & (d.ece_origin == d.host_origin) & d.ece_origin.notna() & (d.ece_origin != "UNMAPPED")
    e_rec = ok.sum() / len(d)
    e_prec = ok.sum() / d["p"].sum() if d["p"].sum() else np.nan

    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.6))

    # A: assembly chimerism vs (misleading) completeness
    bars = ax[0].bar(["chimeric\nE. coli contigs", "CheckM2\ncompleteness"],
                     [chim_pct, comp], color=[RED, GREY], width=0.6)
    for b, v in zip(bars, [chim_pct, comp]):
        ax[0].annotate(f"{v:.0f}%", (b.get_x() + b.get_width() / 2, v),
                       textcoords="offset points", xytext=(0, 4), ha="center", fontsize=11)
    ax[0].set(ylim=(0, 108), ylabel="percent",
              title="A. Strain mixture -> chimeric assembly\n(8 E. coli strains, 96.6-98.6% ANI)")

    # B: de-novo strain-level linkage (recall vs precision)
    bars = ax[1].bar(["recall", "precision"], [e_rec, e_prec], color=[BLUE, ORANGE], width=0.6)
    for b, v in zip(bars, [e_rec, e_prec]):
        ax[1].annotate(f"{v:.2f}", (b.get_x() + b.get_width() / 2, v),
                       textcoords="offset points", xytext=(0, 4), ha="center", fontsize=11)
    ax[1].set(ylim=(0, 1.1), ylabel="strain-level metric",
              title="B. De-novo ECE->strain linkage\nstays precise on chimeric contigs")

    fig.suptitle("MODIFI is robust to strain-mixture assembly chimerism (E. coli panel)",
                 fontsize=12.5, y=1.03)
    fig.tight_layout()
    out = f"{OUT}/ecoli_panel_assembly_linkage.pdf"
    fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150, bbox_inches="tight")
    print(f"wrote {out} | chimeric={chim_pct:.0f}% completeness={comp:.0f}% | "
          f"denovo rec/prec={e_rec:.2f}/{e_prec:.2f}")


if __name__ == "__main__":
    main()
