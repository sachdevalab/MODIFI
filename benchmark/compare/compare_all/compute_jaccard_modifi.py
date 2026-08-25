"""Pairwise Jaccard on DEREPLICATED motif sets, using MODIFI's own MotifFilter
(scripts/derep_motifs.py). Each tool's motifs are pooled across contigs, grouped
by identifier, and dereplicated with MotifFilter; two dereplicated motifs (from
different tools) are matched by MotifFilter's own relatedness (subset / reverse
complement / shared core). Jaccard = shared clusters / total clusters."""
import csv, os, sys
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/scripts")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/compare/compare_all")
from Bio.Seq import Seq
from derep_motifs import MotifFilter
from meta_motif_comparison import OUT, FIG, selected_contigs, MIN_FRAC, MIN_SITES

MF = MotifFilter([])  # for its relatedness/core methods
def rc(s): return str(Seq(s).reverse_complement())

DIRS = {"MODIFI": f"{OUT}/modifi_subread_216/motifs",
        "ipdSummary": f"{OUT}/motifs_compare/percontig/ipdSummary",
        "fibertools": f"{OUT}/motifs_compare/percontig/fibertools"}

def rows_of(path):
    out = []
    if not os.path.isfile(path): return out
    for r in csv.DictReader(open(path)):
        m = (r.get("motifString") or "").strip().upper()
        if not m: continue
        try:
            f=float(r.get("fraction",0) or 0); nd=int(float(r.get("nDetected",0) or 0))
            ng=int(float(r.get("nGenome",0) or 0)); cp=int(float(r.get("centerPos",0) or 0))
        except ValueError: continue
        if f<MIN_FRAC or nd<MIN_SITES: continue
        out.append({"motif":m,"centerPos":cp,"host_meth":nd,"host_total":ng if ng else nd,"indentifier":f"{m}_{cp}"})
    return out

def derep_tool(d):
    allrows=[]
    for c in selected_contigs(): allrows += rows_of(f"{d}/{c}.motifs.csv")
    if not allrows: return []
    df=pd.DataFrame(allrows).groupby("indentifier",as_index=False).agg(
        {"motif":"first","centerPos":"first","host_meth":"sum","host_total":"sum"})
    return MotifFilter(df.to_dict("records")).filter()

def related(a, b):
    m1, m2 = a["motif"], b["motif"]
    if MF.is_subset_or_reverse_complement(m1, m2): return True
    c1, c2 = MF.extract_core_simple(m1), MF.extract_core_simple(m2)
    return c1 == c2 or c1 == rc(c2)

def jaccard(A, B):
    parent = {("A",i):("A",i) for i in range(len(A))}
    parent.update({("B",j):("B",j) for j in range(len(B))})
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for i in range(len(A)):
        for j in range(len(B)):
            if related(A[i], B[j]): parent[find(("A",i))]=find(("B",j))
    comp={}
    for node in parent: comp.setdefault(find(node),[0,0])[0 if node[0]=="A" else 1]+=1
    shared=sum(1 for a,b in comp.values() if a and b)
    return shared/len(comp) if comp else 0.0

reps = {t: derep_tool(d) for t, d in DIRS.items()}
for t in reps: print(f"{t}: {len(reps[t])} dereplicated motifs")
names = list(reps)
with open(f"{FIG}/fig_jaccard_drep.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["tool1","tool2","J_dereplicated"])
    for i in range(len(names)):
        for j in range(len(names)):
            J = 1.0 if i==j else jaccard(reps[names[i]], reps[names[j]])
            w.writerow([names[i], names[j], round(J,4)])
            if i<j: print(f"{names[i]} vs {names[j]}: J={J:.3f}")
