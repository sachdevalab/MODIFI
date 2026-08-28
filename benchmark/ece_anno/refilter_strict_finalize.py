#!/usr/bin/env python3
"""Step 3: apply the revised artifact gate (zero-coverage fraction > 0.10) to the step-1
candidates and write the final revised filter-passing ECE set + counts."""
import csv
from collections import defaultdict

O = "/home/shuaiw/borg/revision/ece_anno/expanded"
RF = f"{O}/refilter_strict"
ZCF_MAX = 0.10   # pass requires zero-coverage fraction <= 10%

# zero-cov lookup
zc = {}
for r in csv.DictReader(open(f"{RF}/zerocov.csv")):
    try: zc[(r["sample"], r["contig"])] = float(r["zero_cov_frac"])
    except: zc[(r["sample"], r["contig"])] = None

step1 = list(csv.DictReader(open(f"{O}/refilter_strict_step1.csv")))
counts = defaultdict(lambda: defaultdict(int))
final = []
n_missing_cov = 0
for r in step1:
    if r["candidate"] != "1":
        continue
    key = (r["sample"], r["MGE"]); typ = r["MGE_type"]
    counts["candidate"][typ] += 1
    z = zc.get(key)
    if z is None:
        n_missing_cov += 1
        flag_art = True          # no coverage evidence -> treat as artifact (conservative)
    else:
        flag_art = z > ZCF_MAX
    r["zero_cov_frac"] = "" if z is None else f"{z:.4f}"
    r["flag_artifact"] = int(flag_art)
    r["pass_revised"] = int(not flag_art)
    if not flag_art:
        counts["final"][typ] += 1
        final.append(r)

# write final passing set
cols = list(step1[0].keys()) + ["zero_cov_frac", "flag_artifact", "pass_revised"]
cols = list(dict.fromkeys(cols))
with open(f"{O}/filterpass_revised.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
    for r in final: w.writerow({k: r.get(k, "") for k in cols})

def line(stage, label):
    p, v = counts[stage]["plasmid"], counts[stage]["virus"]
    print(f"  {label:34} plasmid {p:5}  virus {v:5}  total {p+v:5}")
print("REVISED-CRITERIA final filter-passing set")
line("candidate", "candidates (P2&P3&not-chromosomal)")
print(f"  candidates missing coverage (->artifact): {n_missing_cov}")
line("final", "FINAL PASS (zero-cov frac <= 10%)")
print(f"written: {O}/filterpass_revised.csv ({len(final)} rows)")
