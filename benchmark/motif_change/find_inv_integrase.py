#!/usr/bin/env python3
"""Find proteins in a genomic region (GFF) and extract their sequences from Prokka .faa."""

from pathlib import Path

from Bio import SeqIO

GFF = Path("/home/shuaiw/borg/paper/E_faecalis/prokka/infant_14_31_C.gff")
FAA = Path("/home/shuaiw/borg/paper/E_faecalis/prokka/infant_14_31_C.faa")
REGION_START = 2323857
REGION_END = 2327116


def cds_in_region(gff_path, start, end):
    """Yield (locus_id, feat_start, feat_end) for CDS overlapping [start, end]."""
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            feat_start = int(parts[3])
            feat_end = int(parts[4])
            if feat_end < start or feat_start > end:
                continue
            attrs = dict(p.split("=", 1) for p in parts[8].split(";") if "=" in p)
            locus = attrs.get("ID") or attrs.get("locus_tag")
            if locus:
                yield locus, feat_start, feat_end


if __name__ == "__main__":
    loci = list(cds_in_region(GFF, REGION_START, REGION_END))
    records = SeqIO.to_dict(SeqIO.parse(FAA, "fasta"))

    for locus, s, e in loci:
        if locus in records:
            print(records[locus].format("fasta"), end="")
        else:
            print(f"# {locus} [{s}-{e}] not in .faa", file=__import__("sys").stderr)
