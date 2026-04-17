#!/usr/bin/env python3
"""
Klingon Detection Script

Maps Klingon reference sequences onto assembly contigs using minimap2
and identifies potential Klingon sequences based on alignment quality and coverage.
Results from all samples are merged into a single output file.
"""

import subprocess
import os
import sys
import pandas as pd
import argparse
from pathlib import Path

KLINGON_REF = "/home/shuaiw/borg/paper/natasha/klingon/all.klingo.fasta"
PREFIX_TABLE = "/home/shuaiw/MODIFI/assembly_pipe/prefix_table_soil.tab"
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))


def run_minimap2(assembly_fasta, paf_output, klingon_ref, threads=32):
    for file_path in [assembly_fasta, klingon_ref]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Input file not found: {file_path}")

    cmd = ["minimap2", "-t", str(threads), "-x", "asm5", "--secondary=no",
           klingon_ref, assembly_fasta]

    with open(paf_output, "w") as outfile:
        result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"✗ Minimap2 failed for {assembly_fasta}:\n{result.stderr}")
        return False
    return True


def parse_paf_file(paf_file, min_identity=0.8, min_coverage=0.5):
    columns = [
        "query_name", "query_length", "query_start", "query_end",
        "strand", "target_name", "target_length", "target_start", "target_end",
        "matches", "alignment_length", "mapping_quality",
    ]

    alignments = []
    with open(paf_file) as f:
        for line in f:
            if line.strip():
                fields = line.strip().split("\t")
                if len(fields) >= 12:
                    alignments.append(fields[:12])

    if not alignments:
        return pd.DataFrame()

    df = pd.DataFrame(alignments, columns=columns)
    numeric_cols = ["query_length", "query_start", "query_end",
                    "target_length", "target_start", "target_end",
                    "matches", "alignment_length", "mapping_quality"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col])

    df["identity"] = df["matches"] / df["alignment_length"]
    df["query_coverage"] = (df["query_end"] - df["query_start"]) / df["query_length"]
    df["target_coverage"] = (df["target_end"] - df["target_start"]) / df["target_length"]

    return df[(df["identity"] >= min_identity) & (df["query_coverage"] >= min_coverage)]


def find_assembly():
    fasta_dict = {}
    for line in open(PREFIX_TABLE):
        fields = line.strip().split()
        if len(fields) >= 2:
            fasta_dict[fields[1]] = fields[0]
    return fasta_dict


def main():
    parser = argparse.ArgumentParser(
        description="Find Klingon sequences across all soil assemblies using minimap2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--klingon_ref", default=KLINGON_REF)
    parser.add_argument("--prefix_table", default=PREFIX_TABLE)
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--min_identity", type=float, default=0.8)
    parser.add_argument("--min_coverage", type=float, default=0.5)
    args = parser.parse_args()

    fasta_dict = {}
    for line in open(args.prefix_table):
        fields = line.strip().split()
        if len(fields) >= 2:
            fasta_dict[fields[1]] = fields[0]

    all_summaries = []

    for prefix, assembly_fasta in sorted(fasta_dict.items()):
        print(f"\n🔍 Processing {prefix}...")
        paf_output = os.path.join(OUTPUT_DIR, f"{prefix}.minimap2.paf")

        ok = run_minimap2(assembly_fasta, paf_output, args.klingon_ref, args.threads)
        if not ok:
            continue

        df = parse_paf_file(paf_output, args.min_identity, args.min_coverage)
        if df.empty:
            print(f"  No hits for {prefix}")
            continue

        for target_name, group in df.groupby("target_name"):
            best = group.loc[group["identity"].idxmax()]
            all_summaries.append({
                "sample_name":      prefix,
                "seq_name":         best["query_name"],
                "klingon_ref":      target_name,
                "identity":         round(best["identity"], 3),
                "query_coverage":   round(best["query_coverage"], 3),
                "target_coverage":  round(best["target_coverage"], 3),
                "alignment_length": best["alignment_length"],
                "length":           best["target_length"],
            })
        print(f"  ✓ {len(df.groupby('target_name'))} Klingon hits in {prefix}")

    if not all_summaries:
        print("\n⚠ No Klingon contigs found across any sample.")
        return

    summary_df = pd.DataFrame(all_summaries).sort_values(["sample_name", "identity"], ascending=[True, False])

    summary_file = os.path.join(OUTPUT_DIR, "klingon_contigs_summary.tsv")
    summary_df.to_csv(summary_file, index=False, sep="\t")
    print(f"\n✓ Combined summary saved to: {summary_file}")

    list_file = os.path.join(OUTPUT_DIR, "klingon_contigs.txt")
    with open(list_file, "w") as f:
        f.write("seq_name\n")
        for name in summary_df["seq_name"].unique():
            f.write(f"{name}\n")
    print(f"✓ Klingon reference list saved to: {list_file}")
    print(f"\n🎉 Done. {len(summary_df)} total hits across {summary_df['sample_name'].nunique()} samples.")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
