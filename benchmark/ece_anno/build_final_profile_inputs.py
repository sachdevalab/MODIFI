#!/usr/bin/env python3
"""Prep inputs for the final strict-set (filterpass_FINAL, 3,880) profile figure.
Outputs (under borg/revision/ece_anno/expanded/final_profile/):
  ece_table.csv       : per-ECE sample,MGE,type,environment,length,depth,best_score,best_spec,linked
  linkage_table.csv   : per high-conf linkage w/ RESCUED host: MGE,type,env,host,host_phylum,
                        MGE_gc,host_gc,cos_sim,MGE_cov,host_cov
  crispr_by_sample.csv: Sample,MGE Type,Consistent Linkages  (rescued host, spacer mismatch 0)
"""
import sys, os, re, csv, glob
from collections import defaultdict
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/network")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/isolation")
sys.path.insert(0, "/home/shuaiw/MODIFI/scripts")
import plot_linkage_data as pld
from sample_object import get_detail_taxa_name, classify_taxa
import pysam
from get_kmer_freq import Calc_kmer_freq

A = "/home/shuaiw/borg/revision/ece_anno"; L = "/home/shuaiw/borg/revision/ece_linkage"
RUN2 = "/home/shuaiw/borg/paper/run2"
OUT = f"{A}/expanded/final_profile"; os.makedirs(OUT, exist_ok=True)
ctg_taxa = pld.get_ctg_taxa(RUN2 + "/")
def classified(c): return not re.search("Unclassified", get_detail_taxa_name(ctg_taxa.get(c, "NA")))
def phylum(c): return classify_taxa(ctg_taxa.get(c, "NA"), "phylum")

def num(x):
    try: return float(x)
    except: return None

# final set
fin = {(r["sample"], r["MGE"]): r for r in csv.DictReader(open(f"{A}/expanded/filterpass_FINAL.csv"))}
# depth + environment
depth = {}; env = {}
for r in csv.DictReader(open(f"{A}/expanded/combined_metagenome_eces.csv")):
    depth[(r["sample"], r["MGE"])] = num(r.get("MGE_cov")); env[(r["sample"], r["MGE"])] = r.get("environment")

# per-sample top-host summary row (has GC/cos/cov) + best score lookup
top_row = {}   # key -> host_summary row (top host)
best = {}      # key -> (fs, sp) best across sources
def add_best(s, m, fs, sp):
    fs = num(fs); sp = num(sp)
    if fs is None or sp is None: return
    k = (s, m)
    if k not in best or fs > best[k][0]: best[k] = (fs, sp)
for hs in glob.glob(f"{L}/*/meth_all/host_summary.csv") + glob.glob(f"{L}/ocean_1/meth_all_c*/host_summary.csv") + glob.glob(f"{L}/*/meth_missing/host_summary.csv") + glob.glob(f"{L}/*/meth/host_summary.csv"):
    s = hs.split("/ece_linkage/")[1].split("/")[0]
    if "ocean_1/meth_all_c" in hs: s = "ocean_1"
    for r in csv.DictReader(open(hs)):
        k = (s, r["MGE"]); add_best(s, r["MGE"], r.get("final_score"), r.get("specificity"))
        if k not in top_row: top_row[k] = r   # first (from per_sample this is the top host)
for r in csv.DictReader(open(f"{A}/expanded/genomad_passing_linked_scores_all.csv")):
    add_best(r["sample"], r["seq_name"], r.get("final_score"), r.get("specificity"))
for r in csv.DictReader(open(f"{L}/new_ece_linkage_all.csv")):
    add_best(r["sample"], r["MGE"], r.get("final_score"), r.get("specificity"))

def hp(s, mge):
    for d in ("meth_all", "meth_missing", "meth"):
        p = f"{L}/{s}/{d}/hosts/{mge}.host_prediction.csv"
        if os.path.exists(p): return p
    g = glob.glob(f"{L}/{s}/meth_all_c*/hosts/{mge}.host_prediction.csv")
    return g[0] if g else None

# --- ECE table (all 3,880) ---
with open(f"{OUT}/ece_table.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["sample","MGE","type","environment","length","depth","best_score","best_spec","linked_highconf"])
    for (s, m), r in fin.items():
        fs, sp = best.get((s, m), (0, 1))
        hc = int(fs > 0.8 and sp < 0.001)
        w.writerow([s, m, r["MGE_type"], env.get((s,m),""), r.get("mge_len"), depth.get((s,m),""), fs, sp, hc])

# --- rescued linkage table (high-conf) ---
kc = Calc_kmer_freq()
fa_cache = {}
def fa(s):
    if s not in fa_cache:
        fa_cache[s] = pysam.FastaFile(f"{RUN2}/{s}/{s}.hifiasm.p_ctg.rename.fa")
    return fa_cache[s]
mdep_cache = {}
def mdep(s):
    if s not in mdep_cache:
        d = {}
        for mm in ("methylation4", "methylation3"):
            p = f"{RUN2}/{s}/{s}_{mm}/mean_depth.csv"
            if os.path.exists(p):
                for line in open(p):
                    f = line.strip().split(",");
                    try: d[f[0]] = float(f[1])
                    except: pass
                break
        mdep_cache[s] = d
    return mdep_cache[s]
def gc_of(seq):
    seq = seq.upper(); n = len(seq)
    return (seq.count("G") + seq.count("C")) / n if n else 0
def kmer_norm(seq):
    r = kc.get_kmer_count(seq.upper())
    return r
def cos_pair(s, mge, host):
    try:
        a = fa(s).fetch(mge); b = fa(s).fetch(host)
    except Exception:
        return None, None, None
    ga, gb = gc_of(a), gc_of(b)
    try:
        cs = kc.get_ctg_sim(a.upper(), b.upper())   # get_kmer_count returns a (counts, norm) tuple;
    except Exception:                                # get_ctg_sim does the normalise+cosine correctly
        cs = None
    return ga, gb, cs

rows = []; n_resc = 0
for (s, m), r in fin.items():
    fs, sp = best.get((s, m), (0, 1))
    if not (fs > 0.8 and sp < 0.001): continue
    tr = top_row.get((s, m))
    top_host = tr["host"] if tr else None
    # choose rescued host
    if top_host and classified(top_host):
        host = top_host; src = "top"
    else:
        p = hp(s, m); host = None
        if p:
            cand = []
            for h in csv.DictReader(open(p)):
                f2 = num(h["final_score"]); s2 = num(h["specificity"])
                if f2 is not None and s2 is not None and f2 > 0.8 and s2 < 0.001: cand.append((f2, h["host"]))
            cand.sort(reverse=True)
            host = next((h for f2, h in cand if classified(h)), None)
        if host is None: continue   # no classified host at any rank -> excluded (matches network)
        src = "rescued"; n_resc += 1
    # metrics
    if src == "top" and tr:
        mgc, hgc, cs = num(tr.get("MGE_gc")), num(tr.get("host_gc")), num(tr.get("cos_sim"))
        mcov, hcov = num(tr.get("MGE_cov")), num(tr.get("host_cov"))
    else:
        mgc, hgc, cs = cos_pair(s, m, host)
        md = mdep(s); mcov, hcov = md.get(m), md.get(host)
    rows.append(dict(sample=s, MGE=m, type=r["MGE_type"], environment=env.get((s,m),""),
                     host=host, host_phylum=phylum(host), host_src=src,
                     MGE_gc=mgc, host_gc=hgc, cos_sim=cs, MGE_cov=mcov, host_cov=hcov))
with open(f"{OUT}/linkage_table.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"ECE table: {len(fin)} | linkage table: {len(rows)} (rescued non-top: {n_resc})")

# --- CRISPR consistency (rescued host) ---
def read_spacer(sample):
    f = f"/home/shuaiw/borg/revision/network/spacer/{sample}/{sample}_mge_spacer_hits.filter.tsv"
    d = defaultdict(set)
    if os.path.exists(f):
        for r in csv.DictReader(open(f), delimiter="\t"):
            tgt, qc, fm = r.get("target_id"), r.get("query_contig_id"), r.get("full_mismatch")
            if not tgt or not qc: continue
            try:
                if int(float(fm)) > 0: continue
            except (TypeError, ValueError): pass
            d[tgt].add(qc)
    return d
spc = {}; counts = defaultdict(lambda: defaultdict(int)); tot = 0; val = 0
byt = defaultdict(int); bytv = defaultdict(int)
for r in rows:
    s = r["sample"]
    if s not in spc: spc[s] = read_spacer(s)
    tot += 1; byt[r["type"]] += 1
    if r["MGE"] in spc[s] and r["host"] in spc[s][r["MGE"]]:
        val += 1; bytv[r["type"]] += 1; counts[s][r["type"]] += 1
with open(f"{OUT}/crispr_by_sample.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["Sample", "MGE Type", "Consistent Linkages"])
    for s in sorted(counts):
        for t in ("plasmid", "virus"):
            if counts[s][t]: w.writerow([s, t, counts[s][t]])
print(f"CRISPR: {val}/{tot} consistent ({100*val/tot:.1f}%) | plasmid {bytv['plasmid']}/{byt['plasmid']}, virus {bytv['virus']}/{byt['virus']}")
print(f"outputs -> {OUT}/")
