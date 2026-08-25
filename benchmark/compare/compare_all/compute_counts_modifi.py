"""Count motifs using MODIFI's own dereplication (scripts/derep_motifs.py
MotifFilter: subset / reverse-complement / core-based collapse), the same code
the pipeline ships. Writes per-contig counts (panels a-c) and global
distinct+non-redundant counts (panel d)."""
import csv, os, sys
import pandas as pd
sys.path.insert(0, "/home/shuaiw/MODIFI/scripts")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/compare/compare_all")
from derep_motifs import MotifFilter
from meta_motif_comparison import OUT, FIG, selected_contigs, MIN_FRAC, MIN_SITES

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
            f = float(r.get("fraction", 0) or 0); nd = int(float(r.get("nDetected", 0) or 0))
            ng = int(float(r.get("nGenome", 0) or 0)); cp = int(float(r.get("centerPos", 0) or 0))
        except ValueError: continue
        if f < MIN_FRAC or nd < MIN_SITES: continue
        out.append({"motif": m, "centerPos": cp, "host_meth": nd,
                    "host_total": ng if ng else nd, "indentifier": f"{m}_{cp}"})
    return out

def derep_count(rows):
    if not rows: return 0
    return len(MotifFilter(list(rows)).filter())

cs = selected_contigs()

# ---- per-contig (panels a-c) ----
with open(f"{FIG}/fig_percontig_counts.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["contig", "MODIFI", "ipdSummary", "fibertools"])
    for c in cs:
        w.writerow([c] + [derep_count(rows_of(f"{d}/{c}.motifs.csv")) for d in DIRS.values()])
print("wrote per-contig counts (MotifFilter dereplicated)")

# ---- global (panel d): pool across contigs, group by identifier, then derep ----
with open(f"{FIG}/fig_bars.counts.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["tool", "n_motifs", "n_nonredundant"])
    for t, d in DIRS.items():
        allrows = []
        for c in cs:
            allrows += rows_of(f"{d}/{c}.motifs.csv")
        if allrows:
            df = pd.DataFrame(allrows).groupby("indentifier", as_index=False).agg(
                {"motif": "first", "centerPos": "first", "host_meth": "sum", "host_total": "sum"})
            distinct = len(df)
            pooled = df.to_dict("records")
            nr = len(MotifFilter(pooled).filter())
        else:
            distinct = nr = 0
        w.writerow([t, distinct, nr]); print(f"{t}: distinct={distinct} non-redundant={nr}")
