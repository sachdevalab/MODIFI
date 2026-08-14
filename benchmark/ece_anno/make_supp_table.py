#!/usr/bin/env python3
"""Build the metagenome ECE supplementary table: all 379 ECEs with annotation +
loose-criteria include/exclude decision + a per-row exclusion reason, plus a column
dictionary. Writes an .xlsx (two sheets) if openpyxl is available, else .tsv + a
sidecar column-legend .tsv."""
import os
import pandas as pd

DATA = "/home/shuaiw/borg/revision/ece_anno/metagenome"
CSV = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv"
OUT_DIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
OUT_XLSX = os.path.join(OUT_DIR, "ece_supp_table.xlsx")
OUT_TSV = os.path.join(OUT_DIR, "ece_supp_table.tsv")
OUT_COLS = os.path.join(OUT_DIR, "ece_supp_table_columns.tsv")

# Column dictionary (ordered) -> description
COLUMNS = [
    ("sample", "Metagenome sample the ECE was assembled from"),
    ("seq_name", "ECE contig identifier"),
    ("type", "geNomad element type (plasmid or virus)"),
    ("host", "Linked host contig (from the MODIFI linkage table)"),
    ("host_taxa", "GTDB species-level annotation of the host"),
    ("host_phylum", "GTDB phylum of the host"),
    ("environment", "Habitat of the sample"),
    ("mge_len", "ECE contig length (bp)"),
    # completeness (P1)
    ("circular", "Completeness: closed cycle in the assembly graph (contig ends in C) or a DTR"),
    ("complete_linear", "Completeness: inverted terminal repeat (ITR); complete linear genome"),
    ("genomad_topology", "geNomad topology field (No terminal repeats / DTR / ITR / Provirus)"),
    # geNomad robustness (P2)
    ("genomad_score", "geNomad calibrated classification score"),
    ("genomad_fdr", "geNomad false-discovery rate for the classification"),
    ("survives_strict", "geNomad robustness: score >= 0.8 AND FDR <= 0.05"),
    # plasmid marker (P3)
    ("conjscan_type", "ConjScan call: Conjugative / Mobilizable / dCONJ / None"),
    ("genomad_conjugation_genes", "Count of geNomad-annotated conjugation genes"),
    ("plas_rep_initiation", "Count of plasmid replication-initiation (Rep) proteins (Pfam)"),
    ("plas_partition", "Count of partition proteins ParA/ParB (Pfam)"),
    # viral marker (P3)
    ("vir_terminase_large", "Pfam hits: large terminase subunit"),
    ("vir_terminase_small", "Pfam hits: small terminase subunit"),
    ("vir_major_capsid", "Pfam hits: major capsid protein"),
    ("vir_portal", "Pfam hits: portal protein"),
    ("vir_n_classes", "Distinct structural classes present (terminase L/S, capsid, portal) via Pfam OR VOGdb"),
    ("vog_hallmark_hits", "Total VOGdb hits on the element (informational)"),
    # chromosomal negative
    ("scmg_count", "Distinct bacterial/archaeal single-copy core genes hit (chromosomal signal; used only in strict mode)"),
    ("scmg_fraction", "scmg_count as a fraction of the marker set"),
    ("rrna_count", "rRNA genes (16S/23S) detected on the contig; presence flags chromosomal origin"),
    # coverage
    ("cov_mean", "ECE read depth (copy number numerator)"),
    ("host_depth", "Host read depth"),
    ("cov_ratio", "ECE-to-host depth ratio (used only in strict mode)"),
    ("cov_cv", "Per-base coefficient of variation of ECE coverage (uniformity)"),
    # decision
    ("support_genomad", "Positive line P2 passed (geNomad robustness)"),
    ("support_marker", "Positive line P3 passed (element-type structural/mobility marker)"),
    ("support_circular", "Positive line P1 passed (completeness); recorded, not required in loose mode"),
    ("support_coverage", "Positive line P4 passed (coverage distinctness); recorded, not required in loose mode"),
    ("flag_chromosomal", "Negative flag: likely chromosomal (loose = rRNA present)"),
    ("flag_artifact", "Negative flag: likely assembly artifact (low/erratic depth, no support)"),
    ("retained", "yes if very-high-confidence under the loose metagenome criteria, else no"),
    ("exclusion_reason", "For excluded ECEs, which loose requirement failed"),
]


def exclusion_reason(r):
    if r["very_high_confidence"]:
        return ""
    reasons = []
    if not r["support_genomad"]:
        reasons.append("fails geNomad robustness (score<0.8 or FDR>0.05)")
    if not r["support_marker"]:
        reasons.append("no element-type marker"
                       if r["type"] != "virus" else "<2 viral structural classes")
    if r.get("flag_chromosomal"):
        reasons.append("rRNA gene present")
    if r.get("flag_artifact"):
        reasons.append("assembly-artifact coverage")
    return "; ".join(reasons)


def main():
    e = pd.read_csv(os.path.join(DATA, "ece_evidence_all.tsv"), sep="\t")
    meta = pd.read_csv(CSV)[["MGE", "host", "host_phylum", "environment"]].rename(
        columns={"MGE": "seq_name"}).drop_duplicates("seq_name")
    df = e.merge(meta, on="seq_name", how="left")
    df["retained"] = df["very_high_confidence"].map({True: "yes", False: "no"})
    df["exclusion_reason"] = df.apply(exclusion_reason, axis=1)
    ordered = [c for c, _ in COLUMNS if c in df.columns]
    df = df[["sample"] + [c for c in ordered if c != "sample"]]

    cols_df = pd.DataFrame(COLUMNS, columns=["column", "description"])
    try:
        with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
            df.to_excel(xw, sheet_name="ECEs", index=False)
            cols_df.to_excel(xw, sheet_name="Column_definitions", index=False)
        print(f"wrote {OUT_XLSX} ({len(df)} ECEs, {df['retained'].eq('yes').sum()} retained)")
    except Exception as ex:
        df.to_csv(OUT_TSV, sep="\t", index=False)
        cols_df.to_csv(OUT_COLS, sep="\t", index=False)
        print(f"openpyxl unavailable ({ex}); wrote {OUT_TSV} + {OUT_COLS}")
    # always emit the legend as a standalone tsv too
    cols_df.to_csv(OUT_COLS, sep="\t", index=False)
    print(f"retained {int(df['retained'].eq('yes').sum())}/{len(df)} "
          f"(plasmid {int(((df.type=='plasmid')&(df.retained=='yes')).sum())}, "
          f"virus {int(((df.type=='virus')&(df.retained=='yes')).sum())})")


if __name__ == "__main__":
    main()
