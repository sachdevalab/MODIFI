#!/usr/bin/env python3
"""Build RefSeq-mapping inputs for the K. pneumoniae ECE set, RECOMPUTED from the NEW strict network.
K. pneumoniae set = ECEs whose linkage host is Klebsiella pneumoniae in the strict 317-linkage table,
mapped to the new 95%-ANI MGE clusters; query = all members of those clusters present in all_mge_revised.fa.
Writes under borg/revision/kp_eces/refseq_map/:
  kp_plasmids.fna / kp_viruses.fna  (cluster members split by MGE type)
  kp_members_meta.tsv               (MGE, MGE_cluster, MGE_type, is_kp_linked, inferred_host_genus)"""
import sys, csv
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/network")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/isolation")
import plot_linkage_data as pld
from sample_object import classify_taxa

FP = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile/linkage_table.csv"
MGECLU = "/home/shuaiw/borg/revision/network/MGE_cluster/megablast.cluster.95ani.tsv"
FA = "/home/shuaiw/borg/revision/network/all_mge_revised.fa"
TYPES = "/home/shuaiw/borg/revision/ece_anno/expanded/filterpass_revised_final.csv"
OUT = "/home/shuaiw/borg/revision/kp_eces/refseq_map"

ctg_taxa = pld.get_ctg_taxa("/home/shuaiw/borg/paper/run2/")
mge_clu, clu_mge = pld.read_mge_cluster(MGECLU)
lk = list(csv.DictReader(open(FP)))

# per-contig MGE type from the 3,976-ECE filter table
mtype = {}
for r in csv.DictReader(open(TYPES)):
    mtype[r["MGE"]] = r["MGE_type"]

# K. pneumoniae-linked ECEs -> clusters + inferred genus
kp_linked = set(); cl_type = {}; kp_gen = {}
kp_clusters = set()
for r in lk:
    sp = classify_taxa(ctg_taxa.get(r["host"], "NA"), "species")
    if "Klebsiella pneumoniae" in str(sp):
        c = mge_clu.get(r["MGE"], r["MGE"])
        kp_clusters.add(c); kp_linked.add(r["MGE"])
        cl_type[c] = r["type"]
        g = classify_taxa(ctg_taxa.get(r["host"], "NA"), "genus").replace("g__", "")
        kp_gen[r["MGE"]] = g

# all cluster members present in the ECE fasta
hdr = set()
for l in open(FA):
    if l[0] == ">": hdr.add(l[1:].split()[0])
# a cluster member is counted ONLY if it is a high-confidence linked ECE (in the 317 linkage set),
# not every co-member of the 95%-ANI cluster from the full ECE set.
linked_set = {r["MGE"] for r in lk}
members = []
for c in sorted(kp_clusters):
    for m in clu_mge.get(c, [c]):
        if m in hdr and m in linked_set:
            t = mtype.get(m, cl_type.get(c, "plasmid"))
            members.append((m, c, t))

# write meta + split fasta
mmeta = {m: (c, t) for m, c, t in members}
with open(f"{OUT}/kp_members_meta.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["MGE", "MGE_cluster", "MGE_type", "is_kp_linked", "inferred_host_genus"])
    for m, c, t in members:
        w.writerow([m, c, t, int(m in kp_linked), kp_gen.get(m, "")])

plas = {m for m, c, t in members if t == "plasmid"}
viru = {m for m, c, t in members if t == "virus"}
def write_fa(names, path):
    keep = False; n = 0
    with open(path, "w") as out:
        for line in open(FA):
            if line[0] == ">":
                keep = line[1:].split()[0] in names
                if keep: n += 1
            if keep: out.write(line)
    return n
npl = write_fa(plas, f"{OUT}/kp_plasmids.fna")
nvi = write_fa(viru, f"{OUT}/kp_viruses.fna")
print(f"K. pneumoniae clusters (new network): {len(kp_clusters)}; linked ECEs: {len(kp_linked)}")
print(f"cluster members in fasta: {len(members)} (plasmid {len(plas)}, virus {len(viru)})")
print(f"  wrote kp_plasmids.fna ({npl}), kp_viruses.fna ({nvi}), kp_members_meta.tsv")
print("  clusters:", ", ".join(sorted(kp_clusters)))
