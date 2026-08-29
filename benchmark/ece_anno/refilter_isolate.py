#!/usr/bin/env python3
"""Apply the strict ECE criteria to the FULL geNomad isolate ECE set (4052). Isolate-specific:
P2 = geNomad-strict ONLY, and P1 completeness (circular OR complete_linear) is REQUIRED (isolate
genomes are complete, so circularity/terminal-repeat signals are reliable).

Gates (all must pass):
  P1 completeness : circular OR complete_linear
  P2 geNomad      : genomad_score>=0.8 & genomad_fdr<=0.05
  P3 marker       : virus >=2 grouped {terminase(L/S),major_capsid,portal} (Pfam u VOGdb);
                    plasmid ConjScan!=None OR geNomad conjugation>0 OR plasmid Rep/Par Pfam>0
  not chromosomal : rrna_count==0 & scmg_count<5
  not artifact    : zero_cov_frac<=0.10
  size/depth      : length>=5000 & depth>=5
"""
import csv, glob, os
from collections import defaultdict

A = "/home/shuaiw/borg/revision/ece_anno"
ISO = f"{A}/isolate_all"
DB = f"{A}/db/viral_markers"
GS_MIN, GF_MAX, SCMG_MIN, ZCF_MAX = 0.8, 0.05, 5, 0.10

def f(x):
    try: return float(x)
    except: return None

# --- grouped viral classes per contig (Pfam UNION VOGdb tblouts under isolate_all) ---
acc2class = {r["accession"].split(".")[0]: r["class"]
             for r in csv.DictReader(open(f"{DB}/marker_map.tsv"), delimiter="\t")}
vog2class = {r["vog"]: r["class"]
             for r in csv.DictReader(open(f"{DB}/vog_structural_map.tsv"), delimiter="\t")}
present = defaultdict(set)
g2c = lambda g: g.rsplit("_", 1)[0]
for tbl in glob.glob(f"{ISO}/**/viral_pfam.tblout", recursive=True):
    for ln in open(tbl):
        if ln.startswith("#") or not ln.strip(): continue
        p = ln.split(); c = acc2class.get(p[3].split(".")[0])
        if c: present[g2c(p[0])].add(c)
for tbl in glob.glob(f"{ISO}/**/viral_vogdb.tblout", recursive=True):
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

# --- combined ECE set (all geNomad isolate ECEs) ---
comb = {}
for r in csv.DictReader(open(f"{ISO}/isolate_all_eces.csv")):
    comb[(r["prefix"], r["mge_contig"])] = {"type": r["mge_type"],
        "len": f(r["mge_length"]), "depth": f(r["mge_depth"])}

# --- evidence ---
ev = {}
for pth in glob.glob(f"{ISO}/per_sample/*/ece_evidence.tsv"):
    s = pth.split("/per_sample/")[1].split("/")[0]
    for r in csv.DictReader(open(pth), delimiter="\t"):
        ev[(s, r["seq_name"])] = r

# --- zero-coverage ---
zcf = {}
zp = f"{ISO}/zerocov.csv"
if os.path.exists(zp):
    for r in csv.DictReader(open(zp)):
        zcf[(r["sample"], r["contig"])] = f(r["zero_cov_frac"])

counts = defaultdict(lambda: defaultdict(int)); rows_out = []
for key, c in comb.items():
    s, seq = key; typ = c["type"]; e = ev.get(key)
    # P1 completeness
    circ = (e.get("circular") in ("1", "True", "true")) if e else False
    clin = (e.get("complete_linear") in ("1", "True", "true")) if e else False
    p1 = circ or clin
    # P2 geNomad-strict
    gs = f(e["genomad_score"]) if e else None; gf = f(e["genomad_fdr"]) if e else None
    p2 = gs is not None and gs >= GS_MIN and gf is not None and gf <= GF_MAX
    # P3 marker by type
    if typ == "virus":
        gn = grouped_n(seq); p3 = gn >= 2; marker = f"vir_grouped={gn}"
    else:
        ct = e.get("conjscan_type") if e else None
        conj = ct not in (None, "None", "", "nan")
        gcg = (f(e.get("genomad_conjugation_genes")) or 0) if e else 0
        pmt = (f(e.get("plas_marker_total")) or 0) if e else 0
        p3 = conj or gcg > 0 or pmt > 0; marker = f"conj={int(conj)},gcg={int(gcg)},pmt={int(pmt)}"
    # negatives + size/depth
    rrna = (f(e["rrna_count"]) if e else 0) or 0
    scmg = (f(e["scmg_count"]) if e else 0) or 0
    flag_chr = rrna > 0 or scmg >= SCMG_MIN
    z = zcf.get(key)
    flag_art = (z is not None and z > ZCF_MAX) or (z is None)     # missing cov = conservative artifact
    ln, dp = c["len"], c["depth"]
    ok_size = ln is not None and ln >= 5000 and dp is not None and dp >= 5

    counts["combined"][typ] += 1
    if p1: counts["P1"][typ] += 1
    if p1 and p2: counts["P1_P2"][typ] += 1
    if p1 and p2 and p3: counts["P1_P2_P3"][typ] += 1
    if p1 and p2 and p3 and not flag_chr: counts["notchr"][typ] += 1
    if p1 and p2 and p3 and not flag_chr and not flag_art: counts["notart"][typ] += 1
    passed = p1 and p2 and p3 and not flag_chr and not flag_art and ok_size
    if passed: counts["FINAL"][typ] += 1
    rows_out.append({"sample": s, "MGE": seq, "MGE_type": typ, "mge_len": ln, "mge_depth": dp,
        "circular": int(circ), "complete_linear": int(clin), "p1_complete": int(p1),
        "genomad_score": gs, "genomad_fdr": gf, "p2": int(p2), "p3": int(p3), "p3_detail": marker,
        "rrna_count": rrna, "scmg_count": scmg, "flag_chromosomal": int(flag_chr),
        "zero_cov_frac": z, "flag_artifact": int(flag_art), "pass": int(passed)})

with open(f"{ISO}/filterpass_isolate_FINAL.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)

def line(st, lb):
    p, v = counts[st]["plasmid"], counts[st]["virus"]
    print(f"  {lb:34} plasmid {p:5} virus {v:5} total {p+v:5}")
print("ISOLATE strict filter (geNomad-only P2, +P1 completeness)")
line("combined", "combined geNomad pool")
line("P1", "after P1 completeness")
line("P1_P2", "  & P2 geNomad-strict")
line("P1_P2_P3", "  & P3 element marker")
line("notchr", "  & not chromosomal")
line("notart", "  & not artifact (zerocov<=10%)")
line("FINAL", "  & length>=5kb & depth>=5x  [FINAL]")
print(f"wrote {ISO}/filterpass_isolate_FINAL.csv")
