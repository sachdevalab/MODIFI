#!/usr/bin/env python3
"""Build inputs for the IMG-catalogue host cross-check of the FINAL strict linkage set (317 linked ECEs).
Writes, under borg/revision/ece_anno/expanded/img_map/:
  linked_plasmids.fna / linked_viruses.fna  (the linked ECE contigs, split by MGE type)
  linked_eces_meta.tsv                       (MGE, MGE_type, host_genus)  -- inferred host genus per ECE
Host genus = GTDB genus of the linkage host contig (classify_taxa on run2 ctg taxonomy)."""
import sys, csv
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/network")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/isolation")
import plot_linkage_data as pld
from sample_object import classify_taxa

FP = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile/linkage_table.csv"
FA = "/home/shuaiw/borg/revision/network/all_mge_revised.fa"
OUT = "/home/shuaiw/borg/revision/ece_anno/expanded/img_map"

lk = list(csv.DictReader(open(FP)))
ctg_taxa = pld.get_ctg_taxa("/home/shuaiw/borg/paper/run2/")

def genus(host):
    g = classify_taxa(ctg_taxa.get(host, "NA"), "genus")
    if g in ("Unknown", "g__", "", None) or len(g) <= 3:
        return ""
    return g.replace("g__", "")

# meta + per-type MGE lists
mtype = {}; hg = {}
plas, viru = set(), set()
for r in lk:
    m, t = r["MGE"], r["type"]
    mtype[m] = t; hg[m] = genus(r["host"])
    (plas if t == "plasmid" else viru).add(m)

with open(f"{OUT}/linked_eces_meta.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["MGE", "MGE_type", "host_genus"])
    for m in mtype: w.writerow([m, mtype[m], hg[m]])

# split the fasta by type (single pass; keep only linked contigs)
def write_fasta(names, path):
    keep = False; n = 0
    with open(path, "w") as out:
        for line in open(FA):
            if line[0] == ">":
                keep = line[1:].split()[0] in names
                if keep: n += 1
            if keep: out.write(line)
    return n

npl = write_fasta(plas, f"{OUT}/linked_plasmids.fna")
nvi = write_fasta(viru, f"{OUT}/linked_viruses.fna")
ng = sum(1 for m in hg if hg[m])
print(f"linked ECEs: {len(mtype)} (plasmid {len(plas)}, virus {len(viru)})")
print(f"  wrote linked_plasmids.fna ({npl}), linked_viruses.fna ({nvi})")
print(f"  genus-level inferred host: {ng}/{len(mtype)}  (plasmid {sum(1 for m in plas if hg[m])}, virus {sum(1 for m in viru if hg[m])})")
print(f"  meta -> {OUT}/linked_eces_meta.tsv")
