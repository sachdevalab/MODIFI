#!/usr/bin/env python3
"""
Map all RNA samples in reads_directory.srvp.2024.csv to a reference using:
  bowtie2 -1 reads_1 -2 reads_2 -X 2000 --local -x <index> | shrinksam | sambam > out.shrink.sort.bam

Usage:
  python rna_map.py [--ref REF_FASTA] [--index INDEX_PREFIX] [--out_dir DIR] [--csv CSV] [--dry_run]
  If --index is not set, uses REF_FASTA with .fa/.fasta stripped as the bowtie2 index prefix.
"""

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CSV = SCRIPT_DIR / "reads_directory.srvp.2024.csv"
DEFAULT_REF = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.fasta"
DEFAULT_OUT_DIR = "/home/shuaiw/borg/paper/borg_data/RNA/bam"


def main():
    ap = argparse.ArgumentParser(description="Map RNA reads to ref with bowtie2 | shrinksam | sambam")
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV, help="CSV with sample_id, reads_1, reads_2")
    ap.add_argument("--ref", type=str, default=DEFAULT_REF, help="Reference FASTA (used to derive index prefix if --index not set)")
    ap.add_argument("--index", type=str, default=None, help="Bowtie2 index prefix (-x). Default: ref path with .fa/.fasta stripped")
    ap.add_argument("--out_dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for BAMs")
    ap.add_argument("--threads", "-p", type=int, default=20, help="Bowtie2 threads (default: 20)")
    ap.add_argument("--dry_run", action="store_true", help="Only print commands, do not run")
    args = ap.parse_args()

    index_prefix = args.index
    if index_prefix is None:
        index_prefix = args.ref
        for ext in (".fasta", ".fa", ".fna"):
            if index_prefix.lower().endswith(ext):
                index_prefix = index_prefix[: -len(ext)]
                break

    # Bowtie2 -x expects the index prefix; index files are <prefix>.1.bt2, .2.bt2, etc.
    bt2_check = Path(f"{index_prefix}.1.bt2")
    if not bt2_check.exists():
        print(f"Error: Bowtie 2 index not found at prefix: {index_prefix}", file=sys.stderr)
        print("Build the index first with:", file=sys.stderr)
        print(f"  bowtie2-build {shlex.quote(args.ref)} {shlex.quote(index_prefix)}", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(args.csv)
    if "sample_id" not in df.columns or "reads_1" not in df.columns:
        print("CSV must have columns: sample_id, reads_1, reads_2", file=sys.stderr)
        sys.exit(1)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    skipped = 0
    for _, row in df.iterrows():
        sample_id = str(row["sample_id"]).strip()
        r1 = str(row["reads_1"]).strip() if pd.notna(row.get("reads_1")) else ""
        r2 = str(row.get("reads_2", "")).strip() if pd.notna(row.get("reads_2")) else ""
        if not r1:
            skipped += 1
            continue
        # Paired-end: need both reads
        if not r2:
            print(f"Skip {sample_id}: no reads_2 (single-end not run by this script)", file=sys.stderr)
            skipped += 1
            continue
        safe_id = "".join(c if c.isalnum() or c in "._-" else "_" for c in sample_id)
        out_bam = args.out_dir / f"{safe_id}.shrink.sort.bam"

        cmd = f"bowtie2 -1 {shlex.quote(r1)} -2 {shlex.quote(r2)} -X 2000 --local -p {args.threads} -x {shlex.quote(index_prefix)} | shrinksam | sambam > {shlex.quote(str(out_bam))}"
        print(cmd)
        if not args.dry_run:
            ret = subprocess.run(cmd, shell=True)
            if ret.returncode != 0:
                print(f"Failed: {sample_id}", file=sys.stderr)
                sys.exit(ret.returncode)

    print(f"Done. Skipped {skipped} samples.", file=sys.stderr)


if __name__ == "__main__":
    main()
