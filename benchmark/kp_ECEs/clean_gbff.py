#!/usr/bin/env python3
"""Drop redundant bakta 'gene' features from GenBank files so each locus_tag is unique
(LoVis4u treats gene+CDS pairs sharing a locus_tag as duplicate feature ids). Keeps CDS/tRNA/etc."""
import sys, glob, os
from Bio import SeqIO

src, dst = sys.argv[1], sys.argv[2]
os.makedirs(dst, exist_ok=True)
n = 0
for f in sorted(glob.glob(os.path.join(src, "*.gbff"))):
    recs = list(SeqIO.parse(f, "genbank"))
    for r in recs:
        r.features = [ft for ft in r.features if ft.type != "gene"]
        r.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write(recs, os.path.join(dst, os.path.basename(f)), "genbank")
    n += 1
print(f"cleaned {n} gbff -> {dst}")
