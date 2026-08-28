#!/usr/bin/env python3
"""Authoritative revised ECE filter (2026-08-26), matching the criteria flowchart.

P2 confident call (>=1): geNomad STRICT (score>=0.8 & fdr<=0.05) | VS2 STRICT (max_score>=0.9 &
  hallmark>=1) | VS1 (cat1/2/3) | VIBRANT | PlasX>=0.5.
P3 element marker BY MGE_type (not the evidence's own type, to handle type-conflicted ECEs):
  virus  : >=2 of {terminase (large OR small), major_capsid, portal}  [Pfam v37 UNION VOGdb r98]
  plasmid: ConjScan!=None OR geNomad conjugation>0 OR plasmid Rep/Par Pfam>0
flag_chromosomal: rRNA present OR SCMG (bac71/arc76) core-gene count >= 5      [pass: <5 & no rRNA]
flag_artifact   : zero-coverage fraction > 0.10                                 [pass: <=10%]

Emits candidate contig lists (for zero-cov) + step-1 table. Run refilter_v2_finalize.py after
zero-cov is computed for any new candidates.
"""
import csv, glob, os
from collections import defaultdict

A = "/home/shuaiw/borg/revision/ece_anno"
RUN2 = "/home/shuaiw/borg/paper/run2"
REV = "/home/shuaiw/borg/revision/ece_callers"
DB = f"{A}/db/viral_markers"
COMB = f"{A}/expanded/combined_metagenome_eces.csv"
OUT = f"{A}/expanded"
EV_DIRS = ["expanded_new", "expanded_delta", "expanded_refresh", "metagenome_all"]
GENOMAD_SCORE_MIN, GENOMAD_FDR_MAX = 0.8, 0.05
VS2_SCORE_MIN, VS2_HALLMARK_MIN = 0.9, 1
SCMG_MIN = 5

def f(x):
    try: return float(x)
    except: return None

# --- grouped viral classes per contig from tblout union ---
acc2class = {r["accession"].split(".")[0]: r["class"]
             for r in csv.DictReader(open(f"{DB}/marker_map.tsv"), delimiter="\t")}
vog2class = {r["vog"]: r["class"]
             for r in csv.DictReader(open(f"{DB}/vog_structural_map.tsv"), delimiter="\t")}
present = defaultdict(set)
g2c = lambda g: g.rsplit("_", 1)[0]
for tbl in glob.glob(f"{A}/**/viral_pfam.tblout", recursive=True):
    for ln in open(tbl):
        if ln.startswith("#") or not ln.strip(): continue
        p = ln.split(); c = acc2class.get(p[3].split(".")[0])
        if c: present[g2c(p[0])].add(c)
for tbl in glob.glob(f"{A}/**/viral_vogdb.tblout", recursive=True):
    for ln in open(tbl):
        if ln.startswith("#") or not ln.strip(): continue
        p = ln.split(); c = vog2class.get(p[2])
        if c: present[g2c(p[0])].add(c)
def grouped_n(contig):
    cls = present.get(contig, set()); g = 0
    if "terminase_large" in cls or "terminase_small" in cls: g += 1
    if "major_capsid" in cls: g += 1
    if "portal" in cls: g += 1
    return g

# --- combined set ---
comb = {}
for r in csv.DictReader(open(COMB)):
    comb[(r["sample"], r["MGE"])] = {
        "methods": set(m.strip() for m in (r.get("methods") or "").split(",") if m.strip()),
        "type": r["MGE_type"], "len": r.get("mge_len")}

# --- evidence (first non-empty wins) ---
ev = {}
for d in EV_DIRS:
    for pth in glob.glob(f"{A}/{d}/per_sample/*/ece_evidence.tsv"):
        s = pth.split("/per_sample/")[1].split("/")[0]
        for r in csv.DictReader(open(pth), delimiter="\t"):
            ev.setdefault((s, r["seq_name"]), r)

# --- VS2 max_score & hallmark ---
vs2 = {}
for s in sorted(set(k[0] for k in comb)):
    for p in (f"{REV}/{s}/vs2_ok/final-viral-score.tsv",
              f"{RUN2}/{s}/virsorter2/final-viral-score.tsv"):
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                seq = str(r["seqname"]).split("||")[0]; ms, hm = f(r.get("max_score")), f(r.get("hallmark"))
                k = (s, seq)
                if k not in vs2 or (ms or 0) > (vs2[k][0] or 0): vs2[k] = (ms, hm)
            break

cand_by_sample = defaultdict(list); counts = defaultdict(lambda: defaultdict(int)); rows_out = []
for key, c in comb.items():
    s, seq = key; typ = c["type"]; m = c["methods"]; e = ev.get(key)
    gs = f(e["genomad_score"]) if e else None; gf = f(e["genomad_fdr"]) if e else None
    genomad_strict = (any(x.startswith("genomad") for x in m) and gs is not None
                      and gs >= GENOMAD_SCORE_MIN and gf is not None and gf <= GENOMAD_FDR_MAX)
    ms, hm = vs2.get(key, (None, None))
    vs2_strict = ("virsorter2" in m and ms is not None and ms >= VS2_SCORE_MIN
                  and hm is not None and hm >= VS2_HALLMARK_MIN)
    p2 = genomad_strict or vs2_strict or "virsorter1" in m or "vibrant" in m or "plasx" in m

    if typ == "virus":
        gn = grouped_n(seq); p3 = gn >= 2; marker = f"vir_grouped={gn}"
    else:
        ct = e.get("conjscan_type") if e else None
        conj = ct not in (None, "None", "", "nan")
        gcg = (f(e.get("genomad_conjugation_genes")) or 0) if e else 0
        pmt = (f(e.get("plas_marker_total")) or 0) if e else 0
        p3 = conj or gcg > 0 or pmt > 0; marker = f"conj={int(conj)},gcg={int(gcg)},pmt={int(pmt)}"

    rrna = (f(e["rrna_count"]) if e else 0) or 0
    scmg = (f(e["scmg_count"]) if e else 0) or 0
    flag_chr = rrna > 0 or scmg >= SCMG_MIN

    counts["combined"][typ] += 1
    if p2: counts["P2"][typ] += 1
    if p2 and p3: counts["P2_P3"][typ] += 1
    cand = p2 and p3 and not flag_chr
    if cand:
        counts["P2_P3_notchr"][typ] += 1; cand_by_sample[s].append(seq)
    rows_out.append({"sample": s, "MGE": seq, "MGE_type": typ, "mge_len": c["len"],
        "methods": ",".join(sorted(m)), "genomad_strict": int(genomad_strict),
        "vs2_strict": int(vs2_strict), "p2": int(p2), "p3": int(p3), "p3_detail": marker,
        "rrna_count": rrna, "scmg_count": scmg, "flag_chromosomal": int(flag_chr),
        "candidate": int(cand)})

os.makedirs(f"{OUT}/refilter_v2", exist_ok=True)
for s, seqs in cand_by_sample.items():
    open(f"{OUT}/refilter_v2/{s}.candidates.txt", "w").write("\n".join(sorted(set(seqs))) + "\n")
with open(f"{OUT}/refilter_v2_step1.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)

def line(st, lb):
    p, v = counts[st]["plasmid"], counts[st]["virus"]; print(f"  {lb:36} plasmid {p:5} virus {v:5} total {p+v:5}")
print("AUTHORITATIVE revised filter (P3 by MGE_type; grouped terminase; step 1)")
line("combined", "combined pool")
line("P2", "after P2 strict")
line("P2_P3", "after P2 & P3 (type-correct)")
line("P2_P3_notchr", "after not-chromosomal (rRNA/SCMG>=5)")
print(f"candidates for zero-cov: {sum(len(set(v)) for v in cand_by_sample.values())}")
