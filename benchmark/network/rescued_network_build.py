#!/usr/bin/env python3
"""Build the ECE-host network on the FINAL strict set (rescued hosts), and Figure 5d
(largest-component binary heatmap). ECEs clustered at 95% ANI (existing megablast file),
hosts at 99% ANI (new dRep on rescued hosts). Prints all stats for the results text."""
import sys, csv, re, os
from collections import defaultdict, Counter
import networkx as nx
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/network")
sys.path.insert(0, "/home/shuaiw/MODIFI/benchmark/isolation")
import plot_linkage_data as pld
from sample_object import classify_taxa, get_detail_taxa_name

FP = "/home/shuaiw/borg/revision/ece_anno/expanded/final_profile"
N = "/home/shuaiw/borg/revision/network"
NR = f"{N}/rescued_network"
OUT = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
ctg_taxa = pld.get_ctg_taxa("/home/shuaiw/borg/paper/run2/")
def sp_label(c):
    lin = ctg_taxa.get(c, "NA")
    s = classify_taxa(lin, "species")
    if s not in ("Unknown", "s__", None) and len(s) > 3: return s.replace("s__", "")
    return get_detail_taxa_name(lin)

mge_clu, clu_mge = pld.read_mge_cluster(f"{N}/MGE_cluster/megablast.cluster.95ani.tsv")
_, host_clu, clu_host = pld.read_drep_cluster(f"{NR}/dRep_99_out/data_tables/Cdb.csv")

lk = list(csv.DictReader(open(f"{FP}/linkage_table.csv")))
G = nx.DiGraph()
mtype = {}; hlabel = {}; hcluster_of = {}; hc2host = {}
for r in lk:
    mge, host, t = r["MGE"], r["host"], r["type"]
    mc = mge_clu.get(mge, mge)
    hc = host_clu.get(host, host)               # singleton fallback
    mtype[mc] = t
    hcluster_of[host] = hc
    hc2host[hc] = host
    hlabel[hc] = f"{sp_label(host)} ({hc})"
    G.add_node(mc, kind=t); G.add_node(hc, kind="host")
    G.add_edge(mc, hc)
def rep(mc):  # representative ECE contig for a cluster (alphabetically-first member)
    return sorted(clu_mge.get(mc, [mc]))[0]

hosts = [n for n in G if G.nodes[n]["kind"] == "host"]
plas = [n for n in G if G.nodes[n]["kind"] == "plasmid"]
viru = [n for n in G if G.nodes[n]["kind"] == "virus"]
print(f"NETWORK: nodes={G.number_of_nodes()} (host {len(hosts)}, plasmid {len(plas)}, virus {len(viru)}) edges={G.number_of_edges()}")

# ---- enrich node/edge attributes and export the whole network (matches this 411/270 graph) ----
for m in plas + viru:
    G.nodes[m]["type"] = mtype[m]; G.nodes[m]["label"] = rep(m)
for h in hosts:
    mem = hc2host.get(h, h)
    ph = classify_taxa(ctg_taxa.get(mem, "NA"), "phylum")
    G.nodes[h]["type"] = ph if (ph and len(ph) > 3 and ph not in ("Unknown", "p__")) else "Unclassified"
    G.nodes[h]["label"] = hlabel[h]; G.nodes[h]["species"] = hlabel[h].rsplit(" (", 1)[0]
for u, v in G.edges():
    G.edges[u, v]["weight"] = 1; G.edges[u, v]["type"] = mtype.get(u, "")
ND = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised"
os.makedirs(ND, exist_ok=True)
pld.assign_node_colors_for_gml(G)
nx.write_gml(G, f"{ND}/whole_network2.gml")
try:
    pld.write_filtered_gml_abundant_phylum(G, f"{ND}/whole_network2_abundant_phylum.gml", min_host_per_phylum=5)
except Exception as e:
    print("abundant-phylum gml skipped:", e)
pld.export_network_to_excel(G, f"{ND}/network.xlsx")
pld.export_network_to_json(G, f"{ND}/network.json")
print(f"exported whole-network gml/json/xlsx -> {ND}")

ccs = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
lc = ccs[0]
lch = [n for n in lc if G.nodes[n]["kind"] == "host"]
lcp = [n for n in lc if G.nodes[n]["kind"] == "plasmid"]
lcv = [n for n in lc if G.nodes[n]["kind"] == "virus"]
print(f"components={len(ccs)}  largest: host {len(lch)}, plasmid {len(lcp)}, virus {len(lcv)} (total {len(lc)})")

# top-degree host in LC + multi-host plasmid
UG = G.to_undirected()
deg_host = sorted(((UG.degree(h), h) for h in lch), reverse=True)
print("top host clusters by degree (in largest component):")
for d, h in deg_host[:5]: print(f"   {hlabel[h]}  degree={d}")
mge_deg = [(UG.degree(m), m, mtype[m]) for m in (lcp + lcv)]
print("max plasmid multi-host in LC:", max((d for d, m, t in mge_deg if t == "plasmid"), default=0))
print("max virus  multi-host in LC:", max((d for d, m, t in mge_deg if t == "virus"), default=0))
# strain resolution: host species with >1 cluster in LC
sp_clusters = defaultdict(set)
for h in lch: sp_clusters[sp_label([c for c in clu_host[h]][0] if h in clu_host else h)].add(h)
# simpler: species from label
sp2 = defaultdict(set)
for h in lch:
    lbl = hlabel[h].rsplit(" (", 1)[0]; sp2[lbl].add(h)
print("host species with multiple clusters in LC (strain-resolved):")
for s, cs in sp2.items():
    if len(cs) > 1: print(f"   {s}: {len(cs)} clusters")
# environment of LC nodes
print("LC all infant-gut?" )

# ---- Figure 5d: largest-component binary heatmap ----
host_rows = [h for _, h in sorted(((UG.degree(h), h) for h in lch), reverse=True)]
ece_cols = [m for m in (sorted(lcp, key=lambda x:-UG.degree(x)) + sorted(lcv, key=lambda x:-UG.degree(x)))]
def rep(mc):  # representative ECE contig name for label
    members = clu_mge.get(mc, [mc]); return sorted(members)[0]
M = np.zeros((len(host_rows), len(ece_cols)))
for i, h in enumerate(host_rows):
    for j, m in enumerate(ece_cols):
        if G.has_edge(m, h): M[i, j] = 1
fig = plt.figure(figsize=(13, 8))
gs = fig.add_gridspec(2, 2, width_ratios=[len(ece_cols), 6], height_ratios=[3, len(host_rows)], wspace=0.02, hspace=0.02)
axh = fig.add_subplot(gs[1, 0]); axtop = fig.add_subplot(gs[0, 0], sharex=axh); axright = fig.add_subplot(gs[1, 1], sharey=axh)
axh.imshow(M, aspect="auto", cmap="Greys", vmin=0, vmax=1, interpolation="none")
axh.set_yticks(range(len(host_rows))); axh.set_yticklabels([hlabel[h] for h in host_rows], fontsize=7)
axh.set_xticks(range(len(ece_cols)))
axh.set_xticklabels([rep(m) for m in ece_cols], rotation=90, fontsize=5)
for tick, m in zip(axh.get_xticklabels(), ece_cols):
    tick.set_color("#F8766D" if mtype[m] == "plasmid" else "#00BFC4")
axtop.bar(range(len(ece_cols)), M.sum(0), color="grey"); axtop.set_ylabel("degree", fontsize=8); axtop.tick_params(labelbottom=False); axtop.set_xlim(-0.5, len(ece_cols)-0.5)
axright.barh(range(len(host_rows)), M.sum(1), color="grey"); axright.set_xlabel("degree", fontsize=8); axright.invert_yaxis(); axright.tick_params(labelleft=False); axright.set_ylim(len(host_rows)-0.5, -0.5)
axh.text(-0.16, 1.06, "d", transform=axh.transAxes, fontsize=18, fontweight="bold")
# NOTE: the canonical Fig 5d is the crisp R geom_tile version (plot_component_heatmap_revised.R).
# Write this matplotlib version under a distinct name so re-running this builder never clobbers it.
fig.savefig(f"{OUT}/fig5d_largest_component_matplotlib.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig5d_largest_component_matplotlib.png", bbox_inches="tight", dpi=200)
# source data
with open(f"{OUT}/fig5d_largest_component_sourcedata.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["host_cluster", "ECE_cluster", "MGE_type", "ECE_representative"])
    for i, h in enumerate(host_rows):
        for j, m in enumerate(ece_cols):
            if M[i, j]: w.writerow([hlabel[h], m, mtype[m], rep(m)])
print("wrote fig5d_largest_component.{pdf,png} + sourcedata")
