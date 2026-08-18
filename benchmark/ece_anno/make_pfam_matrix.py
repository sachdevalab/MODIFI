#!/usr/bin/env python3
"""Per-ECE x Pfam-gene presence table: for every ECE, the hit count of each curated
Pfam marker profile (11 viral structural + 10 plasmid Rep/Par), parsed from the stored
viral_pfam.tblout / plasmid_pfam.tblout. Writes an .xlsx (Pfam_presence + Pfam_legend
sheets) to tmp/rev_figs/ece_anno/. Default dataset = metagenome; --dataset isolate for
the isolate set."""
import argparse
import glob
import os
import sys

import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import anno_utils as au

DB = "/home/shuaiw/borg/revision/ece_anno/db"
OUT_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
ROOTS = {
    "metagenome": "/home/shuaiw/borg/revision/ece_anno/metagenome",
    "isolate": "/home/shuaiw/borg/revision/ece_anno",
}


def load_legend():
    frames = []
    for kind, sub in [("viral", "viral_markers"), ("plasmid", "plasmid_markers")]:
        m = pd.read_csv(os.path.join(DB, sub, "marker_map.tsv"), sep="\t")
        m["target"] = kind
        frames.append(m)
    leg = pd.concat(frames, ignore_index=True)
    leg["column"] = leg["pfam_name"] + "|" + leg["accession"]
    return leg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=list(ROOTS), default="metagenome")
    args = ap.parse_args()
    root = ROOTS[args.dataset]
    os.makedirs(OUT_DIR, exist_ok=True)

    leg = load_legend()
    acc2col = dict(zip(leg["accession"], leg["column"]))
    col_order = leg["column"].tolist()  # viral first, then plasmid, in map order

    ece = pd.read_csv(os.path.join(root, "ece_evidence_all.tsv"), sep="\t")
    ece_meta = ece[["sample", "seq_name", "type"]].copy()

    # count hits per (contig, accession) from both Pfam tblouts
    counts = {}  # seq_name -> {column: count}
    for tbl in glob.glob(os.path.join(root, "per_sample", "*",
                                      "*_pfam.tblout")):
        df = au.parse_hmmer_tblout(tbl)
        if df.empty:
            continue
        df["contig"] = df["gene_id"].map(au.gene_to_contig)
        for (contig, acc), sub in df.groupby(["contig", "query_acc"]):
            col = acc2col.get(acc)
            if col is None:
                continue
            counts.setdefault(contig, {}).setdefault(col, 0)
            counts[contig][col] += int(sub["gene_id"].nunique())

    mat = pd.DataFrame(0, index=ece_meta["seq_name"], columns=col_order, dtype=int)
    for contig, d in counts.items():
        if contig in mat.index:
            for col, n in d.items():
                mat.at[contig, col] = n
    mat = mat.reset_index()
    out = ece_meta.merge(mat, on="seq_name", how="left").fillna(0)
    for c in col_order:
        out[c] = out[c].astype(int)
    out["n_viral_pfam_hits"] = out[[c for c in col_order if c in
                                    set(leg.loc[leg.target == "viral", "column"])]].sum(axis=1)
    out["n_plasmid_pfam_hits"] = out[[c for c in col_order if c in
                                      set(leg.loc[leg.target == "plasmid", "column"])]].sum(axis=1)

    legend_out = leg[["column", "accession", "pfam_name", "class", "target"]].rename(
        columns={"class": "structural_class"})
    legend_out = pd.concat([legend_out, pd.DataFrame([
        {"column": "n_viral_pfam_hits", "accession": "", "pfam_name": "",
         "structural_class": "sum of the 11 viral columns", "target": "viral"},
        {"column": "n_plasmid_pfam_hits", "accession": "", "pfam_name": "",
         "structural_class": "sum of the 10 plasmid columns", "target": "plasmid"},
    ])], ignore_index=True)

    xlsx = os.path.join(OUT_DIR, f"ece_pfam_genes_{args.dataset}.xlsx")
    try:
        with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
            out.to_excel(xw, sheet_name="Pfam_presence", index=False)
            legend_out.to_excel(xw, sheet_name="Pfam_legend", index=False)
        print(f"wrote {xlsx} ({len(out)} ECEs x {len(col_order)} Pfam columns)")
    except Exception as ex:
        out.to_csv(xlsx.replace(".xlsx", ".tsv"), sep="\t", index=False)
        legend_out.to_csv(xlsx.replace(".xlsx", "_legend.tsv"), sep="\t", index=False)
        print(f"openpyxl unavailable ({ex}); wrote TSV instead")
    # how many ECEs carry >=1 of each column (quick summary)
    hit_any = (out[col_order] > 0).sum().sort_values(ascending=False)
    print("ECEs with >=1 hit per Pfam:")
    print(hit_any.to_string())


if __name__ == "__main__":
    main()
