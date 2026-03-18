#!/usr/bin/env python3
"""
Find plasmids linked to Gram-Positive / Gram-Negative / Archaea hosts (from count_linkages.R logic),
extract their sequences, and merge into a single FASTA under
/home/shuaiw/borg/paper/{gram_positive|gram_negative|archaea}/.
"""
import argparse
import os
import sys
import pandas as pd
from Bio import SeqIO

# Paths (aligned with count_linkages.R and plot_linkage_data.py)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODIFI_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))
PAPER_FIG_DIR = os.path.join(MODIFI_ROOT, "tmp", "figures", "multi_env_linkage", "network_99")
ANNOTATION_PATHS = [
    os.path.join(SCRIPT_DIR, "..", "specificity", "phylum_gram_archaea_annotation.csv"),
    os.path.join(MODIFI_ROOT, "benchmark", "specificity", "phylum_gram_archaea_annotation.csv"),
    os.path.join(SCRIPT_DIR, "phylum_gram_archaea_annotation.csv"),
    os.path.join(MODIFI_ROOT, "benchmark", "specificity", "phylum_gram_archaea_annotation.csv"),
]
ALL_DIR = "/home/shuaiw/borg/paper/run2"
BORG_PAPER_DIR = "/home/shuaiw/borg/paper"

DOMAIN_SUFFIX = {"Gram-Positive": "gram_positive", "Gram-Negative": "gram_negative", "Archaea": "archaea"}
VALID_DOMAINS = list(DOMAIN_SUFFIX.keys())


def find_annotation_file():
    for path in ANNOTATION_PATHS:
        if os.path.exists(path):
            return path
    raise FileNotFoundError("Cannot find phylum_gram_archaea_annotation.csv")


def load_phylum_class_map(annotation_path):
    df = pd.read_csv(annotation_path)
    return dict(zip(df["phylum"], df["classification"]))


def get_ctg_ref(sample: str, mge_id: str) -> str:
    """Contig path for MGE (same logic as My_contig in sample_object.py)."""
    work_dir = os.path.join(ALL_DIR, sample, f"{sample}_methylation4")
    return os.path.join(work_dir, "contigs", f"{mge_id}.fa")


def get_domain_paths(domain: str):
    """Return (out_dir, base_name, out_fasta) for a domain (Gram-Positive | Gram-Negative | Archaea)."""
    suffix = DOMAIN_SUFFIX[domain]
    base_name = f"{suffix}_linked_plasmids"
    out_dir = os.path.join(BORG_PAPER_DIR, suffix)
    out_fasta = os.path.join(out_dir, f"{base_name}.fasta")
    return out_dir, base_name, out_fasta


def main():
    parser = argparse.ArgumentParser(description="Extract plasmid FASTA for plasmids linked to a host domain.")
    parser.add_argument(
        "--domain",
        choices=VALID_DOMAINS,
        default="Gram-Positive",
        help="Host domain: Gram-Positive, Gram-Negative, or Archaea (default: Gram-Positive)",
    )
    args = parser.parse_args()
    domain = args.domain
    out_dir, base_name, out_fasta = get_domain_paths(domain)

    csv_path = os.path.join(PAPER_FIG_DIR, "mge_host_gc_cov.csv")
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.", file=sys.stderr)
        sys.exit(1)

    annotation_path = find_annotation_file()
    phylum_class_map = load_phylum_class_map(annotation_path)

    gc_df = pd.read_csv(csv_path)
    # Strip p__ from host_phylum and map to Gram-Positive / Gram-Negative / Archaea
    gc_df["host_phylum_clean"] = gc_df["host_phylum"].str.replace("^p__", "", regex=True)
    gc_df["domain_type"] = gc_df["host_phylum_clean"].map(
        lambda p: phylum_class_map.get(p, "Gram-Negative")
    )

    # Plasmids linked to the chosen domain only
    plasmid_dom = gc_df[
        (gc_df["MGE_type"] == "plasmid") & (gc_df["domain_type"] == domain)
    ]
    unique_plasmids = plasmid_dom[["sample", "MGE"]].drop_duplicates()

    os.makedirs(out_dir, exist_ok=True)

    records = []
    seen_ids = set()
    missing = []
    for _, row in unique_plasmids.iterrows():
        sample, mge_id = row["sample"], row["MGE"]
        ctg_ref = get_ctg_ref(sample, mge_id)
        if not os.path.exists(ctg_ref):
            missing.append(ctg_ref)
            continue
        for rec in SeqIO.parse(ctg_ref, "fasta"):
            seq_id = rec.id
            if seq_id in seen_ids:
                seq_id = f"{sample}_{seq_id}"
            seen_ids.add(seq_id)
            rec.id = seq_id
            rec.description = f"sample={sample} host_phylum={domain}"
            records.append(rec)

    if missing:
        print(f"Warning: {len(missing)} contig file(s) not found (first 5):", file=sys.stderr)
        for p in missing[:5]:
            print(f"  {p}", file=sys.stderr)

    SeqIO.write(records, out_fasta, "fasta")
    print(f"Wrote {len(records)} plasmid sequence(s) to {out_fasta}")
    print(f"Total unique plasmid–sample pairs ({domain} linked): {len(unique_plasmids)}")


if __name__ == "__main__":
    main()
