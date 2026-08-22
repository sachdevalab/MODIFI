#!/usr/bin/env python3
"""
Build the family-collapsed motif x contig matrix for the JF8 WGA-control run, mirroring the
soil_1 family matrix (same 16 ground-truth families, same kept contigs), from the jf8_wga
motif_profile.csv. Each WGA motif is mapped to one of the 16 families via the curated member
lists (jf8_soil1_family_collapse.csv), with an IUPAC/revcomp fallback for WGA-specific variants.
Family value per contig = max profile over member motifs.

Writes /home/shuaiw/borg/revision/motif_class/jf8_wga_family_matrix.csv
"""
import re, csv, os
import pandas as pd
from Bio.Seq import Seq

D = "/home/shuaiw/borg/revision/motif_class"
J = f"{D}/jf8_wga"
COLLAPSE = "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification/jf8_soil1_family_collapse.csv"

# 16 families in the fixed order used by the existing figure
fam_order = ["GATC","AATCC","GGAGC","CCATC","CAGGAG","CGMAGG","GATGNAG","CAGNNNNNGGA",
             "CATNNNNNGTG","CAGNNNNNCTG","GACNNNNNRGAC","CCANNNNNNCAT","CAAYNNNNNNGTC",
             "GCACNNNNNNGTT","TAANNNNNNCTTG","GAAGNNNNNNNTCC"]

# member -> family (rep) from the collapse file
member2fam = {}
fam_members = {}
with open(COLLAPSE) as fh:
    for line in fh:
        if line.startswith("#") or line.startswith("family_id"): continue
        parts = [p for p in re.split(r"\s{2,}|,", line.strip()) if p]
        if len(parts) < 4: continue
        rep = parts[1]; members = parts[3].split(";")
        fam_members[rep] = members
        for m in members: member2fam[m] = rep

IU = {'A':'A','C':'C','G':'G','T':'T','R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]',
      'M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
def rx(m): return re.compile('^' + ''.join(IU.get(b,b) for b in m) + '$')
def rc(m): return str(Seq(m).reverse_complement())

def same(a, b):
    """same recognition if equal length and IUPAC-compatible on either strand"""
    if len(a) != len(b): return False
    for x, y in ((a, b), (a, rc(b))):
        if rx(x).match(y) or rx(y).match(x): return True
    return False

def assign_family(motif):
    if motif in member2fam: return member2fam[motif]
    # fallback: IUPAC match to any member of any family
    for rep, members in fam_members.items():
        for mem in members:
            if same(motif, mem): return rep
    # fallback 2: the family rep pattern occurs within the motif (or revcomp) as a core
    for rep in fam_order:
        r = rx(rep).pattern[1:-1]
        for cand in (motif, rc(motif)):
            if re.search(r, ''.join(IU.get(b,b) for b in cand)): return rep
    return None

prof = pd.read_csv(f"{J}/motif_profile.csv", index_col=0)
prof.index = [i.rsplit("_", 1)[0] for i in prof.index]   # strip _centerPos
# aggregate to families (max over members)
fam_rows = {}
unmapped = []
for motif, row in prof.groupby(level=0).max().iterrows():
    fam = assign_family(motif)
    if fam is None:
        unmapped.append(motif); continue
    if fam in fam_rows:
        fam_rows[fam] = pd.concat([fam_rows[fam], row], axis=1).max(axis=1)
    else:
        fam_rows[fam] = row
mat = pd.DataFrame({f: fam_rows[f] for f in fam_order if f in fam_rows}).T
mat = mat.reindex(fam_order).fillna(0.0)

# restrict columns to the kept contigs (same as the existing figure)
ann = pd.read_csv(f"{D}/jf8_heatmap_colann.csv")
keep = [c for c in ann["contig"] if c in mat.columns]
mat = mat[keep]
mat.index.name = "motif"
mat.to_csv(f"{D}/jf8_wga_family_matrix.csv")
print(f"WGA family matrix: {mat.shape[0]} families x {mat.shape[1]} contigs")
print(f"unmapped WGA motifs ({len(unmapped)}): {unmapped}")
print(f"per-family max fraction:")
print(mat.max(axis=1).round(3).to_string())
