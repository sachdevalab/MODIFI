#!/usr/bin/env python3
"""Gather ALL geNomad-derived ECEs across the isolate genomes (batch2_results/*/all_mge.tsv) and join
per-sample mean_depth for depth/length. Emits isolate_all_eces.csv for the evidence engine
(cols prefix, mge_contig, mge_type, mge_length, mge_depth, host_depth, host_lineage)."""
import csv, glob, os

B = "/home/shuaiw/borg/paper/isolation/batch2_results"
OUT = "/home/shuaiw/borg/revision/ece_anno/isolate_all"
os.makedirs(OUT, exist_ok=True)

def depth_map(sample):
    for mv in ("_methylation4", "_methylation2", "_methylation"):
        p = f"{B}/{sample}/{sample}{mv}/mean_depth.csv"
        if os.path.exists(p):
            d = {}
            for r in csv.DictReader(open(p)):
                d[r["contig"]] = (r.get("depth", ""), r.get("length", ""))
            return d
    return {}

rows = []
n_pl = n_vi = 0
for mtsv in sorted(glob.glob(f"{B}/*/all_mge.tsv")):
    sample = mtsv.split("/")[-2]
    dm = depth_map(sample)
    for r in csv.DictReader(open(mtsv), delimiter="\t"):
        if "genomad" not in (r.get("methods") or ""):        # geNomad-derived only
            continue
        typ = r["type"]
        if typ not in ("plasmid", "virus"):
            continue
        seq = r["seq_name"]
        dep, ln = dm.get(seq, ("", r.get("length", "")))
        rows.append({"prefix": sample, "mge_contig": seq, "mge_type": typ,
                     "mge_length": ln or r.get("length", ""), "mge_depth": dep,
                     "host_depth": "", "host_lineage": ""})
        n_pl += typ == "plasmid"; n_vi += typ == "virus"

with open(f"{OUT}/isolate_all_eces.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["prefix", "mge_contig", "mge_type", "mge_length",
                                       "mge_depth", "host_depth", "host_lineage"])
    w.writeheader(); w.writerows(rows)
print(f"isolate geNomad ECEs: {len(rows)} (plasmid {n_pl}, virus {n_vi})")
print(f"samples: {len(set(r['prefix'] for r in rows))}")
print(f"with a depth value: {sum(1 for r in rows if r['mge_depth'])}")
print(f"wrote {OUT}/isolate_all_eces.csv")
