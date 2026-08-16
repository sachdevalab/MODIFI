#!/usr/bin/env python3
"""
Emit the table of selected circular soil genomes for the tool comparison.

Selection rule: circular (`_C`) contigs of the test_100 SR-VP assembly that have
BOTH a subread modification result (test_100) AND an aligned per-contig HiFi bam
in run3/bams, ranked by subread depth. Writes genomes.tsv with per-contig
subread depth, HiFi read count and approximate HiFi depth, and the resolved
reference / HiFi-bam paths, so downstream steps and the manuscript can report
depth transparently (HiFi depth is the limiting factor here).

Run:  python select_circular_genomes.py
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys

PREF = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META"
SAMTOOLS = "/shared/software/bin/samtools"
HIFI_BAM_DIR = "/home/shuaiw/borg/bench/soil_zymo/run3/bams"
CIRC_FASTA_DIR = "/home/shuaiw/borg/contigs/circular"
DEPTH_CSV = "/home/shuaiw/borg/paper/ipdsummary/soil_1/modifi.out/test_100/mean_depth.csv"
FAI = "/home/shuaiw/borg/contigs/test_100.fa.fai"
OUT_TSV = "/home/shuaiw/borg/paper/ipdsummary/compare_all/genomes.tsv"

# the ten selected contig-number suffixes (ranked by subread depth, all have HiFi)
SELECTED = ["658", "327", "645", "479", "396", "514", "504", "517", "350", "316"]


def contig_len(fai_path: str) -> dict[str, int]:
    d: dict[str, int] = {}
    with open(fai_path) as f:
        for line in f:
            p = line.split("\t")
            if len(p) >= 2:
                d[p[0]] = int(p[1])
    return d


def subread_depths(csv_path: str) -> dict[str, float]:
    d: dict[str, float] = {}
    with open(csv_path) as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                d[row["contig"]] = float(row["depth"])
            except (KeyError, ValueError):
                continue
    return d


def hifi_reads(bam: str, contig: str) -> int:
    if not os.path.isfile(bam):
        return 0
    try:
        out = subprocess.run(
            [SAMTOOLS, "idxstats", bam], capture_output=True, text=True, check=True
        ).stdout
    except subprocess.CalledProcessError:
        return 0
    total = 0
    for line in out.splitlines():
        p = line.split("\t")
        if p and p[0] == contig:
            total += int(p[2])
    return total


def main() -> int:
    lens = contig_len(FAI)
    sdep = subread_depths(DEPTH_CSV)
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

    fields = [
        "contig", "length", "subread_depth", "hifi_reads", "hifi_depth_approx",
        "ref_fasta", "hifi_bam",
    ]
    rows = []
    for n in SELECTED:
        ctg = f"{PREF}_{n}_C"
        L = lens.get(ctg, 0)
        bam = os.path.join(HIFI_BAM_DIR, f"{ctg}.bam")
        fa = os.path.join(CIRC_FASTA_DIR, f"{ctg}.fasta")
        nreads = hifi_reads(bam, ctg)
        # approx HiFi depth: reads * mean CCS read length (~11 kb) / contig length
        hdep = round(nreads * 11000.0 / L, 1) if L else 0.0
        rows.append({
            "contig": ctg,
            "length": L,
            "subread_depth": sdep.get(ctg, ""),
            "hifi_reads": nreads,
            "hifi_depth_approx": hdep,
            "ref_fasta": fa if os.path.isfile(fa) else "MISSING",
            "hifi_bam": bam if os.path.isfile(bam) else "MISSING",
        })

    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {OUT_TSV} ({len(rows)} genomes)")
    for r in rows:
        print(f"  {r['contig'].split('META_')[1]:9s} len={r['length']:>9} "
              f"subrDep={r['subread_depth']:>8} hifiReads={r['hifi_reads']:>6} "
              f"hifiDep~{r['hifi_depth_approx']:>5} "
              f"ref={'ok' if r['ref_fasta']!='MISSING' else 'MISSING'} "
              f"hifi={'ok' if r['hifi_bam']!='MISSING' else 'MISSING'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
