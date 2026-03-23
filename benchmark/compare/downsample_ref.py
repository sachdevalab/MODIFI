#!/usr/bin/env python3
"""Build nested reference subsets by cumulative contig length (100–500 Mb)."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

raw_ref = "/home/shuaiw/borg/paper/run2/soil_1/soil_1.hifiasm.p_ctg.rename.fa"
mean_depth_csv = (
    "/home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/mean_depth.csv"
)
out_ref_dir = "/home/shuaiw/borg/paper/ipdsummary/subset_ref/"

MIN_CONTIG_BP = 100_000
MIN_DEPTH = 5.0
# Nested targets: each file is a prefix of the next (same contig order as input).
TARGET_MB = (50, 100, 200, 300, 400, 500)
MB = 1_000_000


def load_depth_by_contig(csv_path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with csv_path.open(newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            cid = row["contig"].strip()
            out[cid] = float(row["depth"])
    return out


def samtools_faidx(fasta: Path) -> None:
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)


def main() -> int:
    ref_path = Path(raw_ref)
    depth_path = Path(mean_depth_csv)
    out_dir = Path(out_ref_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    depth_by_id = load_depth_by_contig(depth_path)

    records: list = []
    for rec in SeqIO.parse(ref_path, "fasta"):
        d = depth_by_id.get(rec.id)
        if d is None or d <= MIN_DEPTH:
            continue
        if len(rec) > MIN_CONTIG_BP:
            records.append(rec)

    if not records:
        print(
            "No contigs pass length and depth filters; nothing written.",
            file=sys.stderr,
        )
        return 1

    end_idx = 0
    cum = 0
    for mb in TARGET_MB:
        need = mb * MB
        while end_idx < len(records) and cum < need:
            cum += len(records[end_idx])
            end_idx += 1
        subset = records[:end_idx]
        total_bp = sum(len(r) for r in subset)
        out_path = out_dir / f"test_{mb}.fa"
        SeqIO.write(subset, out_path, "fasta")
        samtools_faidx(out_path)
        print(
            f"Wrote {out_path} + .fai ({len(subset)} contigs, {total_bp} bp; target >= {need} bp)",
            file=sys.stderr,
        )
        if cum < need:
            print(
                f"Warning: only {cum} bp available; target was {need} bp.",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
