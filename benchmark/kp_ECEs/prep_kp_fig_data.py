#!/usr/bin/env python3
"""Prepare the per-cluster table for the updated Fig a/b (14 K. pneumoniae ECE clusters, new network).
Combines kp_clusters_refseq_summary.csv (size, representative mobility, known, AMR/metal gene lists)
with the network host-genus reach (distinct GTDB host genera each cluster links to). Writes
kp_clusters_fig.csv used by fig_kp_combined.R."""
import sys, csv
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/network")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/isolation")
import pandas as pd
import plot_linkage_data as pld
from sample_object import classify_taxa

FIG = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/kp"
LK = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile/linkage_table.csv"
MGECLU = "/home/shuaiw/borg/revision/network/MGE_cluster/megablast.cluster.95ani.tsv"

c = pd.read_csv(f"{FIG}/kp_clusters_refseq_summary.csv", dtype=str).fillna("")

# network host-genus reach per MGE cluster
ctg_taxa = pld.get_ctg_taxa("/home/shuaiw/borg/paper/run2/")
mge_clu, _ = pld.read_mge_cluster(MGECLU)
reach = {}
for r in csv.DictReader(open(LK)):
    cl = mge_clu.get(r["MGE"], r["MGE"])
    g = classify_taxa(ctg_taxa.get(r["host"], "NA"), "genus").replace("g__", "")
    if g and len(g) > 2:
        reach.setdefault(cl, set()).add(g)

import re
def ncount(s):
    return len([x for x in str(s).split(";") if x])
def mobrank(s):  # dominant mobility among the cluster's linked members
    s = str(s)
    if "conjugative" in s: return "conjugative"
    if re.search(r"(^|;)mobilizable", s): return "mobilizable"
    return "non-mobilizable"

out = pd.DataFrame({
    "cluster": c.cluster,
    "type": c.type,
    "n_members": pd.to_numeric(c.n_members),
    "mobility": c.mobility.map(mobrank),
    "known": (pd.to_numeric(c.n_known) > 0).map({True: 1, False: 0}),
    "host_genera": c.cluster.map(lambda x: len(reach.get(x, set()))),
    "n_amr": c.amr_genes.map(ncount),
    "n_metal": c.stress_genes.map(ncount),
})
out = out.sort_values("n_members", ascending=False)
out.to_csv(f"{FIG}/kp_clusters_fig.csv", index=False)
print(out.to_string(index=False))
print(f"\nwrote {FIG}/kp_clusters_fig.csv ({len(out)} clusters)")
