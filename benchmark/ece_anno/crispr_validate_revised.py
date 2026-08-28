#!/usr/bin/env python3
"""Step 5b: CRISPR-consistency of the revised network linkages, reusing the re-BLASTed spacer hits
(revision/network/spacer/<sample>/<sample>_mge_spacer_hits.filter.tsv). A linkage (MGE, host) is
consistent iff a spacer hitting that MGE lives on the linked host contig (mismatch 0), matching
crispr_validation_breakdown.py. Outputs per-linkage validated table + per-sample x MGE_type counts."""
import csv, os
from collections import defaultdict

REV = "/home/shuaiw/borg/revision/network"
GC = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised/mge_host_gc_cov.csv"
OUT_TSV = f"{REV}/crispr_validation_revised.tsv"
OUT_COUNTS = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised/consistent_linkages_by_sample.csv"

def read_spacer(sample, mismatch_allowed=0):
    f = f"{REV}/spacer/{sample}/{sample}_mge_spacer_hits.filter.tsv"
    d = defaultdict(set)
    if not os.path.exists(f):
        return d
    for r in csv.DictReader(open(f), delimiter="\t"):
        tgt, qc, fm = r.get("target_id"), r.get("query_contig_id"), r.get("full_mismatch")
        if not tgt or not qc:            # skip malformed rows
            continue
        try:
            if int(float(fm)) > mismatch_allowed:
                continue
        except (TypeError, ValueError):
            pass                          # keep row if mismatch unparseable (already filtered at BLAST)
        d[tgt].add(qc)
    return d

rows = list(csv.DictReader(open(GC)))
spacer_cache = {}
out = []
counts = defaultdict(lambda: defaultdict(int))   # sample -> type -> n
total = 0; validated = 0
by_type = defaultdict(int); by_type_val = defaultdict(int)
for r in rows:
    s, mge, host, typ = r["sample"], r["MGE"], r["host"], r["MGE_type"]
    if s not in spacer_cache:
        spacer_cache[s] = read_spacer(s)
    sp = spacer_cache[s]
    v = int(mge in sp and host in sp[mge])
    total += 1; validated += v
    by_type[typ] += 1; by_type_val[typ] += v
    if v:
        counts[s][typ] += 1
    out.append({"sample": s, "environment": r.get("environment", ""), "MGE_type": typ,
                "MGE": mge, "host": host, "validated": v})

with open(OUT_TSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(out)

# per-sample x type counts for the figure (only samples with >=1 consistent)
with open(OUT_COUNTS, "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["Sample", "MGE Type", "Consistent Linkages"])
    for s in sorted(counts):
        for t in ("plasmid", "virus"):
            if counts[s][t]:
                w.writerow([s, t, counts[s][t]])

print(f"linkages: {total} | CRISPR-consistent: {validated} ({100*validated/total:.2f}%)")
for t in ("plasmid", "virus"):
    print(f"  {t}: {by_type_val[t]}/{by_type[t]} consistent")
print(f"samples with >=1 consistent linkage: {len(counts)}")
print(f"wrote {OUT_TSV} and {OUT_COUNTS}")
