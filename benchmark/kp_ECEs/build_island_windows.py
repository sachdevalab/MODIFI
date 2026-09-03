#!/usr/bin/env python3
"""Define per-member visualisation windows (locus:start:end:strand) around the resistance islands of
the infant_15_35_C cluster, for focused LoVis4u synteny:
  - metal island : the copper/silver (pco/sil, +/- mer/ars/ter) island (largest resistance cluster
                   containing a pco or sil gene)
  - amr island   : the cluster containing an antibiotic-resistance gene (integron aadA1-sul1-qacEdelta1
                   + ars, or the aph pair)
Writes metal_windows.txt and amr_windows.txt (space-separated tokens).
"""
import csv, sys
from collections import defaultdict

OUT = "/home/shuaiw/borg/revision/kp_eces/synteny"
PAD = 4000
GAP = 20000

lens = {}
with open(f"{OUT}/members_len.tsv") as fh:
    for line in fh:
        if "\t" not in line:
            continue
        p = line.rstrip("\n").split("\t"); lens[p[0]] = int(p[1])

hits = defaultdict(list)
with open(f"{OUT}/synteny_amrfinder.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["Element type"] in ("AMR", "STRESS"):
            hits[r["Contig id"]].append((int(r["Start"]), int(r["Stop"]), r["Gene symbol"], r["Element type"]))


def clusters(v):
    v = sorted(v); cl = [[v[0]]]
    for x in v[1:]:
        if x[0] - cl[-1][-1][1] <= GAP:
            cl[-1].append(x)
        else:
            cl.append([x])
    return cl


def win(m, c):
    s = max(1, min(a for a, b, *_ in c) - PAD)
    e = min(lens[m], max(b for a, b, *_ in c) + PAD)
    return f"{m}:{s}:{e}:1"


metal, amr = [], []
for m, v in sorted(hits.items()):
    cls = clusters(v)
    # metal island = largest cluster containing a pco/sil gene
    metal_cls = [c for c in cls if any(s.startswith(("pco", "sil")) for a, b, s, t in c)]
    if metal_cls:
        metal.append(win(m, max(metal_cls, key=len)))
    # amr island = each cluster containing an antibiotic AMR gene
    for c in cls:
        if any(t == "AMR" for a, b, s, t in c):
            amr.append(win(m, c))

with open(f"{OUT}/metal_windows.txt", "w") as fh:
    fh.write(" ".join(metal))
with open(f"{OUT}/amr_windows.txt", "w") as fh:
    fh.write(" ".join(amr))
print(f"metal island: {len(metal)} members\n  " + "\n  ".join(metal))
print(f"\namr island: {len(amr)} windows\n  " + "\n  ".join(amr))
