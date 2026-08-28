#!/usr/bin/env python3
"""Re-derive the filter-passing ECE set under the REVISED criteria (2026-08-26):

P2 confident call (any of):
  - geNomad STRICT: genomad_score >= 0.8 AND genomad_fdr <= 0.05   (was: any geNomad call)
  - VirSorter2 STRICT: max_score >= 0.9 AND hallmark >= 1          (was: any VS2 call)
  - VirSorter1 (cat1/2/3, unchanged), VIBRANT (non-fragment, unchanged), PlasX >=0.5 (unchanged)
P3 element marker (unchanged): support_marker (virus vir_n_classes>=2; plasmid ConjScan/geNomad-conj/Rep-Par)
flag_chromosomal (REVISED): rRNA present (rrna_count>0) OR SCMG core-gene count >= 5     [pass needs <5 & no rRNA]
flag_artifact  (REVISED): zero-coverage fraction > 0.10  (i.e. <90% of contig covered by reads)  [computed in step 2]

Step 1 (this pass): apply P2 & P3 & flag_chromosomal from precomputed evidence -> candidate set,
and emit per-sample candidate contig lists for the step-2 samtools-coverage job.
"""
import csv, glob, os, sys
from collections import defaultdict

A = "/home/shuaiw/borg/revision/ece_anno"
RUN2 = "/home/shuaiw/borg/paper/run2"
REV = "/home/shuaiw/borg/revision/ece_callers"
COMB = f"{A}/expanded/combined_metagenome_eces.csv"
OUT = f"{A}/expanded"
EV_DIRS = ["expanded_new", "expanded_delta", "expanded_refresh", "metagenome_all"]

GENOMAD_SCORE_MIN, GENOMAD_FDR_MAX = 0.8, 0.05
VS2_SCORE_MIN, VS2_HALLMARK_MIN = 0.9, 1
SCMG_MIN = 5

def f(x):
    try: return float(x)
    except: return None

# --- 1. combined set: key -> (methods set, type, len) ---
comb = {}
for r in csv.DictReader(open(COMB)):
    key = (r["sample"], r["MGE"])
    methods = set(m.strip() for m in (r.get("methods") or "").split(",") if m.strip())
    comb[key] = {"methods": methods, "type": r["MGE_type"], "len": r.get("mge_len")}

# --- 2. evidence: key -> row (first non-empty wins across dirs) ---
ev = {}
for d in EV_DIRS:
    for pth in glob.glob(f"{A}/{d}/per_sample/*/ece_evidence.tsv"):
        sample = pth.split("/per_sample/")[1].split("/")[0]
        for r in csv.DictReader(open(pth), delimiter="\t"):
            key = (sample, r["seq_name"])
            if key not in ev:
                ev[key] = r

# --- 3. VS2 max_score & hallmark: key -> (max_score, hallmark) ---
vs2 = {}
samples = sorted(set(k[0] for k in comb))
for s in samples:
    for p in (f"{REV}/{s}/vs2_ok/final-viral-score.tsv",
              f"{RUN2}/{s}/virsorter2/final-viral-score.tsv"):
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                seq = str(r["seqname"]).split("||")[0]
                ms, hm = f(r.get("max_score")), f(r.get("hallmark"))
                key = (s, seq)
                if key not in vs2 or (ms or 0) > (vs2[key][0] or 0):
                    vs2[key] = (ms, hm)
            break

# --- 4. apply gates ---
cand_by_sample = defaultdict(list)
counts = defaultdict(lambda: defaultdict(int))   # stage -> type -> n
rows_out = []
missing_ev = 0
for key, c in comb.items():
    s, seq = key
    typ = c["type"]; methods = c["methods"]
    e = ev.get(key)
    if e is None:
        missing_ev += 1
    gs = f(e["genomad_score"]) if e else None
    gf = f(e["genomad_fdr"]) if e else None
    genomad_strict = (any(m.startswith("genomad") for m in methods)
                      and gs is not None and gs >= GENOMAD_SCORE_MIN
                      and gf is not None and gf <= GENOMAD_FDR_MAX)
    ms, hm = vs2.get(key, (None, None))
    vs2_strict = ("virsorter2" in methods and ms is not None and ms >= VS2_SCORE_MIN
                  and hm is not None and hm >= VS2_HALLMARK_MIN)
    vs1 = "virsorter1" in methods
    vibrant = "vibrant" in methods
    plasx = "plasx" in methods
    p2 = genomad_strict or vs2_strict or vs1 or vibrant or plasx

    support_marker = (str(e.get("support_marker")).strip().lower() in ("true", "1")) if e else False
    p3 = support_marker

    rrna = f(e["rrna_count"]) if e else 0
    scmg = f(e["scmg_count"]) if e else 0
    flag_chr = (rrna or 0) > 0 or (scmg or 0) >= SCMG_MIN

    counts["combined"][typ] += 1
    if p2: counts["P2"][typ] += 1
    if p2 and p3: counts["P2_P3"][typ] += 1
    cand = p2 and p3 and not flag_chr
    if cand:
        counts["P2_P3_notchr"][typ] += 1
        cand_by_sample[s].append(seq)
    rows_out.append({"sample": s, "MGE": seq, "MGE_type": typ, "mge_len": c["len"],
                     "methods": ",".join(sorted(methods)),
                     "genomad_strict": int(genomad_strict), "vs2_strict": int(vs2_strict),
                     "vs1": int(vs1), "vibrant": int(vibrant), "plasx": int(plasx),
                     "p2": int(p2), "support_marker": int(p3),
                     "rrna_count": rrna or 0, "scmg_count": scmg or 0,
                     "flag_chromosomal": int(flag_chr), "candidate": int(cand)})

# --- write candidate contig lists + full annotation table ---
os.makedirs(f"{OUT}/refilter_strict", exist_ok=True)
for s, seqs in cand_by_sample.items():
    with open(f"{OUT}/refilter_strict/{s}.candidates.txt", "w") as fh:
        fh.write("\n".join(sorted(set(seqs))) + "\n")
with open(f"{OUT}/refilter_strict_step1.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys())); w.writeheader()
    for r in rows_out: w.writerow(r)

# --- report ---
def line(stage, label):
    p, v = counts[stage]["plasmid"], counts[stage]["virus"]
    print(f"  {label:28} plasmid {p:5}  virus {v:5}  total {p+v:5}")
print("REVISED-CRITERIA re-filter (step 1: P2 strict + P3 + chromosomal; artifact pending step 2)")
print(f"combined ECEs: {len(comb)} | evidence rows: {len(ev)} | combined w/o evidence: {missing_ev}")
line("combined", "combined (candidate pool)")
line("P2", "after P2 (strict conf call)")
line("P2_P3", "after P2 & P3 (marker)")
line("P2_P3_notchr", "after not-chromosomal (rRNA/SCMG>=5)")
ncand = sum(len(set(v)) for v in cand_by_sample.values())
print(f"candidates for step-2 coverage: {ncand} across {len(cand_by_sample)} samples")
