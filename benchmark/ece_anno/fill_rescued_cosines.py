#!/usr/bin/env python3
"""Fill the missing tetranucleotide cosine (cos_sim) for the rescued high-confidence linkages, then
propagate the values into the derived per-linkage tables. cos_sim = kc.get_ctg_sim(ECE seq, host seq)."""
import csv, os, sys
import pysam
sys.path.insert(0, "/home/shuaiw/MODIFI/scripts")
from get_kmer_freq import Calc_kmer_freq

RUN2 = "/home/shuaiw/borg/paper/run2"
ECE_FA = "/home/shuaiw/borg/revision/network/all_mge_revised.fa"
LK = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile/linkage_table.csv"
DERIVED = [  # (path, mge_col, sample_col)
    "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/mge_host_gc_cov_final.csv",
    "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised/mge_host_gc_cov.csv",
    "/home/shuaiw/borg/revision/ece_anno/high_conf_linkage/high_conf_linkage_table.csv",
    "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/high_conf_linkage_table.csv",
]
kc = Calc_kmer_freq()
ece_fh = pysam.FastaFile(ECE_FA)
host_fh = {}
def host_seq(sample, host):
    if sample not in host_fh:
        p = f"{RUN2}/{sample}/{sample}.hifiasm.p_ctg.rename.fa"
        host_fh[sample] = pysam.FastaFile(p) if os.path.exists(p) else None
    h = host_fh[sample]
    return h.fetch(host) if (h and host in h.references) else None

def is_blank(v):
    if v in ("", None): return True
    try: return float(v) == 0
    except: return True

# --- fill linkage_table + build a (sample,MGE)->cos map ---
rows = list(csv.DictReader(open(LK)))
cosmap = {}
n = 0
for r in rows:
    if is_blank(r.get("cos_sim")):
        a = ece_fh.fetch(r["MGE"]) if r["MGE"] in ece_fh.references else None
        b = host_seq(r["sample"], r["host"])
        if a and b:
            cs = kc.get_ctg_sim(a.upper(), b.upper())
            r["cos_sim"] = f"{cs:.4f}"; cosmap[(r["sample"], r["MGE"])] = r["cos_sim"]; n += 1
        else:
            print("  MISSING seq for", r["sample"], r["MGE"], r["host"])
with open(LK, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"filled {n} rescued cosines in linkage_table.csv")

# --- propagate into derived tables (update cos_sim where blank, join on sample+MGE) ---
for path in DERIVED:
    if not os.path.exists(path):
        print("  skip (absent):", path); continue
    d = list(csv.DictReader(open(path)))
    if "cos_sim" not in d[0]:
        print("  skip (no cos_sim):", path); continue
    upd = 0
    for r in d:
        k = (r.get("sample"), r.get("MGE"))
        if is_blank(r.get("cos_sim")) and k in cosmap:
            r["cos_sim"] = cosmap[k]; upd += 1
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(d[0].keys())); w.writeheader(); w.writerows(d)
    print(f"  updated {upd} cos_sim in {os.path.basename(path)}")
