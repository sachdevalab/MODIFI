#!/usr/bin/env python3
"""score_e2e.py — score the toy end-to-end run: ECE->host linkage vs assembly completeness.

Joins: MODIFI host_summary (best host contig per ECE) + skani contig_origin (true source
isolate per contig) + CheckM2 per-source-genome completeness + geNomad mge_file + manifest.
An ECE link is CORRECT iff its best host contig traces to the SAME source isolate as the
ECE contig. Reports contig-level recall/precision and linkage success stratified by the
host genome's real assembly completeness — the reviewer's MAG-quality gradient.

Usage: score_e2e.py <C4_toy_dir>
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S0, P0 = 0.5, 0.01   # default operating point


def main():
    D = sys.argv[1]
    hs = pd.read_csv(os.path.join(D, "modifi", "toy", "host_summary.csv"))
    org = pd.read_csv(os.path.join(D, "contig_origin.tsv"), sep="\t")
    o = org.set_index("contig")["origin_sample"].to_dict()
    man = pd.read_csv(os.path.join(D, "toy.manifest.csv"))
    dp = man.set_index("sample")["target_dp"].to_dict()
    sp = man.set_index("sample")["species"].to_dict()
    cm = pd.read_csv(os.path.join(D, "checkm2", "quality_report.tsv"), sep="\t")
    comp = cm.set_index("Name")["Completeness"].to_dict()
    cont = cm.set_index("Name")["Contamination"].to_dict()

    d = hs.copy()
    d["ece_origin"] = d["MGE"].map(o)
    d["host_origin"] = d["host"].map(o)
    d["assigned"] = (d["final_score"] > S0) & (d["specificity"] < P0)
    d["correct"] = d["assigned"] & (d["ece_origin"] == d["host_origin"]) & \
                   d["ece_origin"].notna() & (d["ece_origin"] != "UNMAPPED")
    # completeness / depth of the ECE's TRUE host genome (what the reviewer cares about)
    d["host_completeness"] = d["ece_origin"].map(comp)
    d["host_contamination"] = d["ece_origin"].map(cont)
    d["host_depth"] = d["ece_origin"].map(dp)
    d["host_species"] = d["ece_origin"].map(sp)

    total = len(d)
    mapped = d["ece_origin"].notna() & (d["ece_origin"] != "UNMAPPED")
    n_assigned = int(d["assigned"].sum())
    n_correct = int(d["correct"].sum())
    recall = n_correct / total if total else np.nan
    precision = n_correct / n_assigned if n_assigned else np.nan
    d.to_csv(os.path.join(D, "e2e_predictions.csv"), index=False)

    print(f"\n===== toy end-to-end scoring =====")
    print(f"ECEs (geNomad calls with a MODIFI row): {total}   (origin-mapped: {int(mapped.sum())})")
    print(f"assigned @default(score>{S0},spec<{P0}): {n_assigned}   correct: {n_correct}")
    print(f"contig-level recall={recall:.3f}  precision={precision:.3f}")
    print(f"\nper-ECE (ece -> host, true origins, host completeness):")
    show = ["MGE", "ece_origin", "host", "host_origin", "final_score", "specificity",
            "host_completeness", "host_depth", "correct"]
    print(d[[c for c in show if c in d]].to_string(index=False))

    # stratify linkage success vs the true host genome's completeness
    s = d[mapped].copy()
    if len(s):
        s["comp_bin"] = pd.cut(s["host_completeness"], [0, 50, 70, 90, 100.01],
                               labels=["<50", "50-70", "70-90", "90-100"])
        agg = s.groupby("comp_bin", observed=True).agg(
            n=("correct", "size"), correct=("correct", "sum"),
            mean_depth=("host_depth", "mean")).reset_index()
        agg["success"] = agg["correct"] / agg["n"]
        print(f"\nlinkage success vs host-genome completeness:")
        print(agg.to_string(index=False))

        # figure: success + host completeness/depth
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
        ax[0].scatter(s["host_completeness"], s["final_score"],
                      c=s["correct"].map({True: "#0072B2", False: "#D55E00"}), s=60,
                      edgecolor="k", lw=0.4)
        ax[0].axhline(S0, ls="--", color="#888", lw=1)
        ax[0].set(xlabel="host genome CheckM2 completeness (%)", ylabel="ECE linkage score",
                  title="linkage score vs host completeness\n(blue=correct, orange=miss/wrong)")
        b = agg.dropna(subset=["comp_bin"])
        ax[1].bar(b["comp_bin"].astype(str), b["success"], color="#0072B2")
        ax[1].set(xlabel="host completeness bin (%)", ylabel="fraction ECEs correctly linked",
                  ylim=(0, 1.05), title="linkage success vs completeness")
        fig.tight_layout()
        out = os.path.join(D, "e2e_linkage_vs_completeness.pdf")
        fig.savefig(out, bbox_inches="tight"); fig.savefig(out.replace(".pdf", ".png"), dpi=150)
        print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
