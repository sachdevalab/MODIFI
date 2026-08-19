#!/usr/bin/env python3
"""
Convert SMAC's aggregate 6mA penetrance (per-position, mapped to ref216) into
per-contig modified-base GFFs (with 41-mer context), so the same per-contig
motifMaker step used for the other tools can be applied to SMAC.

SMAC penetrance columns (from 10-Penetrance.pl): contig, pos(1-based), strand,
n_6mA, coverage, penetration(=n_6mA/coverage).

Usage: parse_smac_penetrance.py <penetrance_file> [ref.fa] [out_gff_dir]
  Default ref = ref216.fa, default out dir = motifs_compare/percontig/SMAC_gff/.
  For the per-contig array run, pass the single-contig ref.fa and the contig's
  working dir so it writes just that one <contig>.gff.
Writes: <out_gff_dir>/<contig>.gff
"""
from __future__ import annotations
import os, sys

OUT = "/home/shuaiw/borg/paper/ipdsummary/compare_all_meta"
REF216 = sys.argv[2] if len(sys.argv) > 2 else f"{OUT}/ref216.fa"
GFFDIR = sys.argv[3] if len(sys.argv) > 3 else f"{OUT}/motifs_compare/percontig/SMAC_gff"
FT_FRAC, MINCOV = 0.50, 4      # same thresholds as the other 6mA/aggregate tools

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def _rc(s): return s.translate(_COMP)[::-1]
def _context41(seq, pos1, strand):
    i = pos1 - 1
    lo, hi = i - 20, i + 21
    ctx = "N" * max(0, -lo) + seq[max(0, lo):min(len(seq), hi)] + "N" * max(0, hi - len(seq))
    return _rc(ctx) if strand == "-" else ctx

def load_ref():
    seqs, name, buf = {}, None, []
    for line in open(REF216):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf).upper()
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf).upper()
    return seqs

def main():
    pen = sys.argv[1]
    ref = load_ref()
    os.makedirs(GFFDIR, exist_ok=True)
    handles, n = {}, 0
    for line in open(pen):
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        c, pos, strand = p[0], p[1], p[2]
        try:
            n6, cov, frac = int(p[3]), int(p[4]), float(p[5])
        except ValueError:
            continue
        if c not in ref or cov < MINCOV or frac < FT_FRAC:
            continue
        pos1 = int(pos)
        if strand not in ("+", "-"):
            base = ref[c][pos1 - 1] if 0 <= pos1 - 1 < len(ref[c]) else "N"
            strand = "+" if base == "A" else "-" if base == "T" else "+"
        if c not in handles:
            handles[c] = open(f"{GFFDIR}/{c}.gff", "w")
        ctx = _context41(ref[c], pos1, strand)
        handles[c].write(f"{c}\tSMAC\tmodified_base\t{pos1}\t{pos1}\t{round(frac*100,1)}\t{strand}\t.\t"
                         f"coverage={cov};context={ctx};IPDRatio=2.0\n")
        n += 1
    for h in handles.values():
        h.close()
    print(f"Wrote {len(handles)} per-contig GFFs, {n} 6mA sites (penetration>={FT_FRAC}, cov>={MINCOV}) to {GFFDIR}")

if __name__ == "__main__":
    main()
