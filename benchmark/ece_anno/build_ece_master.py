#!/usr/bin/env python3
"""Assemble one per-ECE master table for the revised filter-passing set (3,976):
sample, seq_name, type, environment, length, depth, final_score, sites, mod_density_per_kb, grp.
length/depth/env from combined_metagenome_eces.csv; final_score from the current linkage lookup
(meth_all + meth_missing + fallback aggregations); mod_density reused from the kept sourcedata
where present, else computed from the ECE's <contig>.gff (modified_base lines, score>=30)."""
import csv, glob, os
from collections import defaultdict

A = "/home/shuaiw/borg/revision/ece_anno"
L = "/home/shuaiw/borg/revision/ece_linkage"
RUN2 = "/home/shuaiw/borg/paper/run2"
MD = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/moddensity_zero_vs_pos_box_sourcedata.csv"
OUT = f"{A}/expanded/ece_master_revised.csv"
MIN_SCORE = 30

def num(x):
    try: return float(x)
    except: return None

# passing set
passing = {(r["sample"], r["MGE"]): r["MGE_type"]
           for r in csv.DictReader(open(f"{A}/expanded/filterpass_revised_final.csv"))}

# combined: length, depth, environment
comb = {}
for r in csv.DictReader(open(f"{A}/expanded/combined_metagenome_eces.csv")):
    comb[(r["sample"], r["MGE"])] = (num(r.get("mge_len")), num(r.get("MGE_cov")), r.get("environment"))

# current linkage score lookup (best final_score per MGE)
score = {}
def add(s, seq, fs, sp):
    fs = num(fs)
    if fs is None: return
    k = (s, seq)
    if k not in score or fs > score[k]: score[k] = fs
for hs in glob.glob(f"{L}/*/meth_all/host_summary.csv"):
    s = hs.split("/ece_linkage/")[1].split("/")[0]
    for r in csv.DictReader(open(hs)): add(s, r["MGE"], r.get("final_score"), None)
for hs in glob.glob(f"{L}/ocean_1/meth_all_c*/host_summary.csv"):
    for r in csv.DictReader(open(hs)): add("ocean_1", r["MGE"], r.get("final_score"), None)
for hs in glob.glob(f"{L}/*/meth_missing/host_summary.csv"):
    s = hs.split("/ece_linkage/")[1].split("/")[0]
    for r in csv.DictReader(open(hs)): add(s, r["MGE"], r.get("final_score"), None)
for r in csv.DictReader(open(f"{A}/expanded/genomad_passing_linked_scores_all.csv")):
    if (r["sample"], r["seq_name"]) not in score: add(r["sample"], r["seq_name"], r.get("final_score"), None)
for r in csv.DictReader(open(f"{L}/new_ece_linkage_all.csv")):
    if (r["sample"], r["MGE"]) not in score: add(r["sample"], r["MGE"], r.get("final_score"), None)

# mod_density reuse
md = {}
for r in csv.DictReader(open(MD)):
    md[(r["sample"], r["seq_name"])] = (num(r.get("sites")), num(r.get("mod_density")))

# gff mod count for missing
def gff_dir(s):
    for mm in ("methylation4", "methylation3"):
        d = f"{RUN2}/{s}/{s}_{mm}/gffs"
        if os.path.isdir(d): return d
    return None
def gff_mod_count(path):
    n = 0
    try:
        with open(path) as fh:
            for ln in fh:
                if ln.startswith("#"): continue
                f = ln.split("\t")
                if len(f) < 6 or f[2] != "modified_base": continue
                try:
                    if float(f[5]) >= MIN_SCORE: n += 1
                except: pass
    except FileNotFoundError:
        return None
    return n

rows = []; n_computed = 0; n_nogff = 0
for (s, mge), typ in passing.items():
    length, depth, env = comb.get((s, mge), (None, None, None))
    fs = score.get((s, mge))
    sites, moddens = md.get((s, mge), (None, None))
    if moddens is None:
        gd = gff_dir(s)
        cnt = gff_mod_count(os.path.join(gd, f"{mge}.gff")) if gd else None
        if cnt is not None and length:
            sites = cnt; moddens = cnt / length * 1000.0; n_computed += 1
        else:
            n_nogff += 1
    grp = None if fs is None else ("score>0" if fs > 0 else "score=0")
    rows.append({"sample": s, "seq_name": mge, "type": typ, "environment": env,
                 "length": length, "depth": depth, "final_score": fs,
                 "sites": sites, "mod_density_per_kb": moddens, "grp": grp})

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"master table: {len(rows)} ECEs -> {OUT}")
print(f"  mod_density reused: {len(rows)-n_computed-n_nogff}, computed: {n_computed}, no gff: {n_nogff}")
print(f"  with final_score: {sum(1 for r in rows if r['final_score'] is not None)}")
print(f"  score=0: {sum(1 for r in rows if r['grp']=='score=0')}, score>0: {sum(1 for r in rows if r['grp']=='score>0')}")
