"""Final dereplicated motif set x tool presence/absence table for manual checking.

Each tool's motifs are pooled across contigs and dereplicated with MODIFI's own
MotifFilter (subset / reverse-complement / shared-core collapse). All three tools'
dereplicated motifs are then clustered together (same MotifFilter relatedness) so
each row is one non-redundant motif; columns mark which tools detected it, and the
*_variants columns list the exact strings each tool contributed (auditable by hand).

Output: tmp/rev_figs/compare_all_meta/motif_presence_by_tool.tsv (next to fig_motif_combined.png).
"""
import csv, os, sys
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/scripts")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/compare/compare_all")
from Bio.Seq import Seq
from derep_motifs import MotifFilter
from meta_motif_comparison import OUT, FIG, selected_contigs, MIN_FRAC, MIN_SITES

MF = MotifFilter([])
def rc(s): return str(Seq(s).reverse_complement())

DIRS = {"MODIFI": f"{OUT}/modifi_subread_216/motifs",
        "ipdSummary": f"{OUT}/motifs_compare/percontig/ipdSummary",
        "fibertools": f"{OUT}/motifs_compare/percontig/fibertools"}
TOOLS = list(DIRS)

def rows_of(path):
    out = []
    if not os.path.isfile(path): return out
    for r in csv.DictReader(open(path)):
        m = (r.get("motifString") or "").strip().upper()
        if not m: continue
        try:
            f = float(r.get("fraction", 0) or 0); nd = int(float(r.get("nDetected", 0) or 0))
            ng = int(float(r.get("nGenome", 0) or 0)); cp = int(float(r.get("centerPos", 0) or 0))
        except ValueError: continue
        if f < MIN_FRAC or nd < MIN_SITES: continue
        out.append({"motif": m, "centerPos": cp, "host_meth": nd,
                    "host_total": ng if ng else nd, "indentifier": f"{m}_{cp}"})
    return out

def derep_tool(d):
    allrows = []
    for c in selected_contigs():
        allrows += rows_of(f"{d}/{c}.motifs.csv")
    if not allrows: return []
    df = pd.DataFrame(allrows).groupby("indentifier", as_index=False).agg(
        {"motif": "first", "centerPos": "first", "host_meth": "sum", "host_total": "sum"})
    return MotifFilter(df.to_dict("records")).filter()

def related(a, b):
    m1, m2 = a["motif"], b["motif"]
    if MF.is_subset_or_reverse_complement(m1, m2): return True
    c1, c2 = MF.extract_core_simple(m1), MF.extract_core_simple(m2)
    return c1 == c2 or c1 == rc(c2)

# 1) per-tool dereplicated reps
reps = {t: derep_tool(d) for t, d in DIRS.items()}
for t in TOOLS:
    sys.stderr.write(f"{t}: {len(reps[t])} dereplicated motifs\n")

# 2) pool as tagged nodes; cluster across all tools by MotifFilter relatedness
nodes = [(t, r) for t in TOOLS for r in reps[t]]
n = len(nodes)
parent = list(range(n))
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for i in range(n):
    for j in range(i + 1, n):
        if related(nodes[i][1], nodes[j][1]):
            parent[find(i)] = find(j)

comp = {}
for i in range(n):
    comp.setdefault(find(i), []).append(i)

# 3) one row per cluster
rowlist = []
for members in comp.values():
    present = {t: 0 for t in TOOLS}
    variants = {t: [] for t in TOOLS}
    best = None
    for idx in members:
        t, r = nodes[idx]
        present[t] = 1
        variants[t].append(r["motif"])
        if best is None or r["host_meth"] > best["host_meth"]:
            best = r
    row = {"motif": best["motif"], "centerPos": best["centerPos"],
           "n_tools": sum(present.values())}
    for t in TOOLS:
        row[t] = present[t]
    for t in TOOLS:
        row[f"{t}_variants"] = ";".join(sorted(set(variants[t])))
    rowlist.append(row)

cols = ["motif", "centerPos", "n_tools"] + TOOLS + [f"{t}_variants" for t in TOOLS]
df = pd.DataFrame(rowlist, columns=cols).sort_values(
    ["n_tools", "motif"], ascending=[False, True]).reset_index(drop=True)

out_tsv = f"{FIG}/motif_presence_by_tool.tsv"
os.makedirs(FIG, exist_ok=True)
df.to_csv(out_tsv, sep="\t", index=False)

print(f"rows (final dereplicated motifs): {len(df)}")
for t in TOOLS:
    print(f"  {t}: {int(df[t].sum())} detected")
print("shared by all 3:", int((df['n_tools'] == 3).sum()))
print("wrote", out_tsv)
