#!/usr/bin/env python3
"""Recompute the viral structural-class count with large+small terminase COLLAPSED into one
'terminase' category (per the intended criteria flowchart): virus P3 = >= 2 of
{terminase, major_capsid, portal}, from Pfam v37 UNION VOGdb r98.

We gather each contig's present-class set from the UNION of all available viral_pfam.tblout +
viral_vogdb.tblout across every evidence dir, so a virus is dropped ONLY if even all evidence
yields < 2 grouped categories. Applies to the revised passing set -> filterpass_revised2.csv.
"""
import csv, glob, os, re
from collections import defaultdict

A = "/home/shuaiw/borg/revision/ece_anno"
DB = f"{A}/db/viral_markers"
PASS = f"{A}/expanded/filterpass_revised.csv"
OUT = f"{A}/expanded/filterpass_revised2.csv"

# class maps
acc2class = {}
for r in csv.DictReader(open(f"{DB}/marker_map.tsv"), delimiter="\t"):
    acc2class[r["accession"].split(".")[0]] = r["class"]
vog2class = {}
for r in csv.DictReader(open(f"{DB}/vog_structural_map.tsv"), delimiter="\t"):
    vog2class[r["vog"]] = r["class"]

def group(cls):
    g = set()
    if "terminase_large" in cls or "terminase_small" in cls: g.add("terminase")
    if "major_capsid" in cls: g.add("major_capsid")
    if "portal" in cls: g.add("portal")
    return g

# contig -> set(raw classes), from union of all tblouts (all dirs)
present = defaultdict(set)
def gene_to_contig(g): return g.rsplit("_", 1)[0]

for tbl in glob.glob(f"{A}/**/viral_pfam.tblout", recursive=True):
    for ln in open(tbl):
        if ln.startswith("#") or not ln.strip(): continue
        f = ln.split()
        gene, acc = f[0], f[3].split(".")[0]
        c = acc2class.get(acc)
        if c: present[gene_to_contig(gene)].add(c)
for tbl in glob.glob(f"{A}/**/viral_vogdb.tblout", recursive=True):
    for ln in open(tbl):
        if ln.startswith("#") or not ln.strip(): continue
        f = ln.split()
        gene, vog = f[0], f[2]
        c = vog2class.get(vog)
        if c: present[gene_to_contig(gene)].add(c)

# apply to passing viruses
rows = list(csv.DictReader(open(PASS)))
kept, dropped = [], []
for r in rows:
    if r["MGE_type"] != "virus":
        r["vir_n_classes_grouped"] = ""
        kept.append(r); continue
    g = group(present.get(r["MGE"], set()))
    r["vir_n_classes_grouped"] = len(g)
    r["grouped_categories"] = "|".join(sorted(g))
    if len(g) >= 2:
        kept.append(r)
    else:
        dropped.append(r)

with open(OUT, "w", newline="") as fh:
    cols = list(rows[0].keys()) + ["vir_n_classes_grouped", "grouped_categories"]
    cols = list(dict.fromkeys(cols))
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader()
    for r in kept: w.writerow(r)

np = sum(1 for r in kept if r["MGE_type"] == "plasmid")
nv = sum(1 for r in kept if r["MGE_type"] == "virus")
print("Grouped-terminase virus rule (>=2 of terminase/capsid/portal):")
print(f"  viruses dropped (grouped <2): {len(dropped)}")
for d in dropped[:60]:
    print(f"    {d['sample']:22} {d['MGE']:28} grouped={d['vir_n_classes_grouped']} [{d.get('grouped_categories','')}]")
print(f"  FINAL revised2: plasmid {np}  virus {nv}  total {np+nv}")
print(f"  written {OUT}")
