#!/usr/bin/env python3
"""
Synthesize a single-contig, motifMaker-ready GFF from fibertools or jasmine
per-position calls: strand inferred from the reference base, plus the 41-mer
`context=` field motifMaker requires (RC on the minus strand). Same logic as the
v1 motif_comparison build, applied per contig.

Usage:
  python pertool_gff.py fibertools <ref.fa> <pileup.bed>      <out.gff> [frac=0.5] [mincov=4]
  python pertool_gff.py jasmine    <ref.fa> <cpg.bed(.gz)>    <out.gff> [pct=50]  [mincov=4]
"""
from __future__ import annotations
import gzip, sys

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def _rc(s): return s.translate(_COMP)[::-1]

def read_one_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()

def ctx41(seq, pos1, strand):
    i = pos1 - 1
    lo, hi = i - 20, i + 21
    c = "N" * max(0, -lo) + seq[max(0, lo):min(len(seq), hi)] + "N" * max(0, hi - len(seq))
    return _rc(c) if strand == "-" else c

def gff_line(name, src, pos1, score, strand, cov, seq):
    return (f"{name}\t{src}\tmodified_base\t{pos1}\t{pos1}\t{score}\t{strand}\t.\t"
            f"coverage={cov};context={ctx41(seq, pos1, strand)};IPDRatio=2.0\n")

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def main():
    tool, ref, calls, out = sys.argv[1:5]
    seq = read_one_fasta(ref)
    name = None
    with open(ref) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]; break
    n = 0
    with open(out, "w") as o, opener(calls) as f:
        if tool == "fibertools":
            frac = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
            mincov = int(sys.argv[6]) if len(sys.argv) > 6 else 4
            for line in f:
                if line[0] == "#":
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9:
                    continue
                cov = int(p[3]); m6a = int(p[8])
                if cov < mincov or cov == 0 or m6a / cov < frac:
                    continue
                pos1 = int(p[1]) + 1
                base = seq[pos1 - 1] if 0 <= pos1 - 1 < len(seq) else "N"
                strand = "+" if base == "A" else "-" if base == "T" else None
                if strand:
                    o.write(gff_line(name, "fibertools", pos1, round(m6a / cov * 100, 1), strand, cov, seq)); n += 1
        elif tool == "jasmine":
            pct = float(sys.argv[5]) if len(sys.argv) > 5 else 50.0
            mincov = int(sys.argv[6]) if len(sys.argv) > 6 else 4
            for line in f:
                if line[0] == "#":
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[4] != "Total":
                    continue
                try:
                    meth = float(p[3]); cov = int(p[5])
                except ValueError:
                    continue
                if cov < mincov or meth < pct:
                    continue
                pos1 = int(p[1]) + 1
                base = seq[pos1 - 1] if 0 <= pos1 - 1 < len(seq) else "N"
                strand = "+" if base == "C" else "-" if base == "G" else None
                if strand:
                    o.write(gff_line(name, "jasmine", pos1, round(meth, 1), strand, cov, seq)); n += 1
        else:
            sys.exit("tool must be fibertools|jasmine")
    print(f"{tool}: wrote {n} sites -> {out}")

if __name__ == "__main__":
    main()
