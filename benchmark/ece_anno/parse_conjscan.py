#!/usr/bin/env python3
"""
Parse CONJScan best_solution_summary.tsv and classify each contig as:
  - conjugative  (has T4SS = complete conjugation machinery)
  - mobilizable  (has MOB relaxase only, no T4SS)
  - dCONJ        (decayed conjugation system, partial T4SS)

Output columns:
  contig, conjugative, mobilizable, dCONJ, conjscan_anno

Optionally parse all_systems.tsv to produce a detailed annotation with:
  gene_name, hit_id, i_eval, score, profile_cov, seq_cov

Usage:
  python parse_conjscan.py <best_solution_summary.tsv> <output.tsv> [all_systems.tsv] [detailed_output.tsv]
"""

import pandas as pd
import sys
import re


def parse_conjscan(input_file, output_file):
    # Skip comment lines starting with #
    df = pd.read_csv(input_file, sep="\t", comment="#")

    # Identify column groups
    data_cols = [c for c in df.columns if c != "replicon"]
    mob_cols  = [c for c in data_cols if "/MOB" in c]
    t4ss_cols = [c for c in data_cols if "/T4SS_" in c]
    dconj_cols = [c for c in data_cols if "/dCONJ_" in c]

    out_rows = []
    for _, row in df.iterrows():
        contig = row["replicon"]

        # Collect hit details per category
        t4ss_hits = []
        for c in t4ss_cols:
            if row[c] > 0:
                # Extract context (Chromosome/Plasmids) and type
                m = re.search(r'CONJScan/(\w+)/T4SS_(\w+)', c)
                if m:
                    t4ss_hits.append(f"T4SS_{m.group(2)}({m.group(1)})")

        mob_hits = []
        for c in mob_cols:
            if row[c] > 0:
                m = re.search(r'CONJScan/(\w+)/MOB', c)
                if m:
                    mob_hits.append(f"MOB({m.group(1)})")

        dconj_hits = []
        for c in dconj_cols:
            if row[c] > 0:
                m = re.search(r'CONJScan/(\w+)/dCONJ_(\w+)', c)
                if m:
                    dconj_hits.append(f"dCONJ_{m.group(2)}({m.group(1)})")

        # Classification
        is_conjugative = len(t4ss_hits) > 0
        is_mobilizable = len(mob_hits) > 0 and not is_conjugative
        is_dconj = len(dconj_hits) > 0

        # Build annotation string
        anno_parts = t4ss_hits + mob_hits + dconj_hits
        conjscan_anno = "; ".join(anno_parts) if anno_parts else ""

        out_rows.append({
            "contig": contig,
            "conjugative": int(is_conjugative),
            "mobilizable": int(is_mobilizable),
            "dCONJ": int(is_dconj),
            "conjscan_anno": conjscan_anno,
        })

    out_df = pd.DataFrame(out_rows)

    # Summary
    n_conj = out_df["conjugative"].sum()
    n_mob  = out_df["mobilizable"].sum()
    n_dconj = out_df["dCONJ"].sum()
    n_any  = ((out_df["conjugative"] + out_df["mobilizable"] + out_df["dCONJ"]) > 0).sum()
    print(f"Total contigs:    {len(out_df)}")
    print(f"Conjugative:      {n_conj}")
    print(f"Mobilizable:      {n_mob}")
    print(f"dCONJ:            {n_dconj}")
    print(f"Any CONJScan hit: {n_any}")

    out_df.to_csv(output_file, sep="\t", index=False)
    print(f"\nOutput saved to: {output_file}")
    return out_df


def parse_detailed_anno(all_systems_file, detailed_output_file):
    """Parse all_systems.tsv to extract per-hit details."""
    rows = []
    with open(all_systems_file) as f:
        for line in f:
            if line.startswith("#") or line.startswith("replicon") or line.strip() == "":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 18:
                continue
            rows.append({
                "gene_name": fields[2],
                "hit_id": fields[1],
                "i_eval": fields[14],
                "score": fields[15],
                "profile_cov": fields[16],
                "seq_cov": fields[17],
            })

    detail_df = pd.DataFrame(rows, columns=["gene_name", "hit_id", "i_eval", "score", "profile_cov", "seq_cov"])
    detail_df = detail_df.sort_values(["hit_id", "gene_name"]).drop_duplicates()

    detail_df.to_csv(detailed_output_file, sep="\t", index=False)
    print(f"Detailed hits: {len(detail_df)}")
    print(f"Detailed output saved to: {detailed_output_file}")
    return detail_df


if __name__ == "__main__":
    if len(sys.argv) not in (3, 5):
        print(f"Usage: {sys.argv[0]} <best_solution_summary.tsv> <output.tsv> [all_systems.tsv] [detailed_output.tsv]")
        sys.exit(1)
    parse_conjscan(sys.argv[1], sys.argv[2])
    if len(sys.argv) == 5:
        parse_detailed_anno(sys.argv[3], sys.argv[4])