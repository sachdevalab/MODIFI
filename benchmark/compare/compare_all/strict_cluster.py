"""Stricter motif clustering + recomputed counts, pairwise Jaccard, and presence table.

Filter: fraction >= 0.4 AND nDetected >= 100.
Two motifs are the SAME only if they have equal length AND one is an IUPAC subset of the
other, OR the reverse complement of one is a subset of the other. No core matching, no
left-aligned different-length matching, no fully-IUPAC-compatible rule. Clusters are the
connected components (single-linkage) under this strict relation.
"""
import csv, os, sys
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/compare/compare_all")
from Bio.Seq import Seq
from meta_motif_comparison import OUT, FIG, selected_contigs

MIN_FRAC, MIN_SITES = 0.4, 100   # stricter site threshold (was 30)

IUPAC = {'A':set('A'),'C':set('C'),'G':set('G'),'T':set('T'),'R':set('AG'),'Y':set('CT'),
         'S':set('GC'),'W':set('AT'),'K':set('GT'),'M':set('AC'),'B':set('CGT'),'D':set('AGT'),
         'H':set('ACT'),'V':set('ACG'),'N':set('ACGT')}
def rc(s): return str(Seq(s).reverse_complement())
def subset(x, y):  # every position: bases(x) subset of bases(y); requires equal length
    return all(IUPAC.get(a, set()) <= IUPAC.get(b, set()) for a, b in zip(x, y))
def strict_related(a, b):
    if len(a) != len(b):
        return False
    if subset(a, b) or subset(b, a):
        return True
    ra = rc(a)
    return subset(ra, b) or subset(b, ra)

DIRS = {"MODIFI": f"{OUT}/modifi_subread_216/motifs",
        "ipdSummary": f"{OUT}/motifs_compare/percontig/ipdSummary",
        "fibertools": f"{OUT}/motifs_compare/percontig/fibertools"}
TOOLS = list(DIRS)

def pool_tool(d):
    agg = {}
    for c in selected_contigs():
        p = f"{d}/{c}.motifs.csv"
        if not os.path.isfile(p): continue
        for r in csv.DictReader(open(p)):
            m = (r.get("motifString") or "").strip().upper()
            if not m: continue
            try:
                f = float(r.get("fraction", 0) or 0); nd = int(float(r.get("nDetected", 0) or 0))
                ng = int(float(r.get("nGenome", 0) or 0)); cp = int(float(r.get("centerPos", 0) or 0))
            except ValueError: continue
            if f < MIN_FRAC or nd < MIN_SITES: continue
            k = (m, cp)
            e = agg.setdefault(k, {"motif": m, "centerPos": cp, "host_meth": 0})
            e["host_meth"] += nd
    return list(agg.values())

def components(nodes):
    """nodes: list of (tool, motifdict). Returns list of member-index lists (single-linkage)."""
    n = len(nodes)
    parent = list(range(n))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    # bucket by length so we only compare equal-length motifs (huge speedup + matches rule)
    from collections import defaultdict
    by_len = defaultdict(list)
    for i, (_, r) in enumerate(nodes):
        by_len[len(r["motif"])].append(i)
    for idxs in by_len.values():
        for ii in range(len(idxs)):
            for jj in range(ii + 1, len(idxs)):
                i, j = idxs[ii], idxs[jj]
                if strict_related(nodes[i][1]["motif"], nodes[j][1]["motif"]):
                    parent[find(i)] = find(j)
    comp = {}
    for i in range(n):
        comp.setdefault(find(i), []).append(i)
    return list(comp.values())

pooled = {t: pool_tool(d) for t, d in DIRS.items()}
for t in TOOLS:
    sys.stderr.write(f"{t}: pooled {len(pooled[t])} motifs (frac>=0.4 & nDetected>=100)\n")

# 1) per-tool non-redundant counts
counts = {t: len(components([(t, r) for r in pooled[t]])) for t in TOOLS}

# 2) pairwise Jaccard (cluster the two tools' pooled motifs together)
def jaccard(A, B):
    nodes = [("A", r) for r in pooled[A]] + [("B", r) for r in pooled[B]]
    comps = components(nodes)
    both = sum(1 for m in comps if any(nodes[i][0] == "A" for i in m) and any(nodes[i][0] == "B" for i in m))
    return both / len(comps) if comps else 0.0
pairs = [("MODIFI", "ipdSummary"), ("MODIFI", "fibertools"), ("ipdSummary", "fibertools")]
jac = {p: jaccard(*p) for p in pairs}

# 3) presence table (all three together)
allnodes = [(t, r) for t in TOOLS for r in pooled[t]]
rows = []
for m in components(allnodes):
    present = {t: 0 for t in TOOLS}; variants = {t: [] for t in TOOLS}; best = None
    for i in m:
        t, r = allnodes[i]; present[t] = 1; variants[t].append(r["motif"])
        if best is None or r["host_meth"] > best["host_meth"]: best = r
    row = {"motif": best["motif"], "centerPos": best["centerPos"], "n_tools": sum(present.values())}
    for t in TOOLS: row[t] = present[t]
    for t in TOOLS: row[f"{t}_variants"] = ";".join(sorted(set(variants[t])))
    rows.append(row)
rows.sort(key=lambda r: (-r["n_tools"], r["motif"]))

cols = ["motif", "centerPos", "n_tools"] + TOOLS + [f"{t}_variants" for t in TOOLS]
with open(f"{FIG}/motif_presence_by_tool.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(f"{FIG}/counts_strict.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["tool", "n_nonredundant"]); [w.writerow([t, counts[t]]) for t in TOOLS]
with open(f"{FIG}/jaccard_strict.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["tool1", "tool2", "J_strict"]); [w.writerow([a, b, round(jac[(a, b)], 4)]) for a, b in pairs]

# ---- figure inputs (strict), in the names fig_combined.R reads ----
# panel d: distinct (pooled) + non-redundant (strict clusters)
with open(f"{FIG}/fig_bars.counts.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["tool", "n_motifs", "n_nonredundant"])
    for t in TOOLS:
        w.writerow([t, len(pooled[t]), counts[t]])
# panel e: pairwise Jaccard (strict) under the column name fig_combined.R expects
with open(f"{FIG}/fig_jaccard_drep.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["tool1", "tool2", "J_dereplicated"])
    for a, b in pairs:
        w.writerow([a, b, round(jac[(a, b)], 4)])
# panels a-c: per-contig non-redundant counts (strict clustering, same filter)
def contig_count(d, c):
    p = f"{d}/{c}.motifs.csv"
    if not os.path.isfile(p): return 0
    agg = {}
    for r in csv.DictReader(open(p)):
        m = (r.get("motifString") or "").strip().upper()
        if not m: continue
        try:
            fr = float(r.get("fraction", 0) or 0); nd = int(float(r.get("nDetected", 0) or 0))
            cp = int(float(r.get("centerPos", 0) or 0))
        except ValueError: continue
        if fr < MIN_FRAC or nd < MIN_SITES: continue
        agg[(m, cp)] = {"motif": m}
    return len(components([("x", v) for v in agg.values()]))
with open(f"{FIG}/fig_percontig_counts.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["contig", "MODIFI", "ipdSummary", "fibertools"])
    for c in selected_contigs():
        w.writerow([c] + [contig_count(DIRS[t], c) for t in TOOLS])
print("wrote fig_bars.counts.csv, fig_jaccard_drep.csv, fig_percontig_counts.csv (strict) to", FIG)

print("\n=== STRICT clustering (equal-length subset or reverse-complement; nDetected>=100) ===")
print("non-redundant motifs:", {t: counts[t] for t in TOOLS})
print("pairwise Jaccard:", {f"{a} vs {b}": round(jac[(a, b)], 3) for a, b in pairs})
print("presence-table rows:", len(rows), "| shared by all 3:", sum(1 for r in rows if r["n_tools"] == 3))
print("wrote motif_presence_by_tool.strict.tsv, counts_strict.csv, jaccard_strict.csv to", FIG)
