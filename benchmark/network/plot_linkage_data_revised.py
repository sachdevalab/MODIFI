#!/usr/bin/env python3
"""Step 4: rebuild the metagenome MGE-host network on the REVISED high-confidence linkage set.

Reuses all helpers from plot_linkage_data.py. Differs only in inputs:
  - per-sample host_summary + all_mge.tsv come from borg/revision/network/per_sample/<sample>/
    (overriding sample_obj.host_sum_file / .mge_file); reference fasta + GTDB stay in run2.
  - MGE clusters  = borg/revision/network/MGE_cluster/megablast.cluster.95ani.tsv (ECE code)
  - host clusters = borg/revision/network/dRep_99_out/data_tables/Cdb.csv (dRep, different code)
Outputs to tmp/rev_figs/ece_anno/network_99_revised/. Non-destructive.
"""
import os, sys
import networkx as nx
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "isolation"))
import plot_linkage_data as pld
from sample_object import My_sample

REV = "/home/shuaiw/borg/revision/network"
RUN2 = "/home/shuaiw/borg/paper/run2/"
META = "/home/shuaiw/MODIFI/assembly_pipe/prefix_table.tab"
MGE_CLU = f"{REV}/MGE_cluster/megablast.cluster.95ani.tsv"
DREP_CLU = f"{REV}/dRep_99_out/data_tables/Cdb.csv"
OUTDIR = "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/network_99_revised/"
os.makedirs(OUTDIR, exist_ok=True)

sample_env = pld.read_metadata(META)
drep_clu_dict, host_clu_dict, clu_host_dict = pld.read_drep_cluster(DREP_CLU)
mge_clu_dict, clu_mge_dict = pld.read_mge_cluster(MGE_CLU)
ctg_taxa_dict = pld.get_ctg_taxa(RUN2)

whole_G = nx.DiGraph()
gc_data = []
cluster_anno_dict = {}
all_host_clu_lineage_dict = {}
host_ctg_set = set()

for prefix in sorted(os.listdir(RUN2)):
    sdir = f"{REV}/per_sample/{prefix}"
    if not os.path.isdir(sdir):
        continue
    if prefix in ["ERR5621427_sludge", "ERR5621429_sludge", "ERR5621430_sludge"]:
        continue
    sample_obj = My_sample(prefix, RUN2)
    # override linkage + MGE-type inputs with the revised per-sample files
    sample_obj.host_sum_file = f"{sdir}/host_summary.csv"
    sample_obj.mge_file = f"{sdir}/all_mge.tsv"
    MGE_type_dict = sample_obj.read_MGE()
    if not MGE_type_dict:
        continue
    G, gc_data, cluster_anno_dict, host_clu_lineage_dict = pld.get_edge(
        cluster_anno_dict, MGE_type_dict, gc_data, sample_env.get(prefix, "NA"),
        prefix, host_clu_dict, mge_clu_dict, sample_obj, ctg_taxa_dict)
    try:
        host_ctg_set = pld.collect_host_ctgs(prefix, sample_obj, ctg_taxa_dict, host_ctg_set)
    except Exception as e:
        print(f"  collect_host_ctgs skipped for {prefix}: {e}")
    whole_G = nx.compose(whole_G, G)
    all_host_clu_lineage_dict.update(host_clu_lineage_dict)

print(f"NETWORK: nodes={whole_G.number_of_nodes()} edges={whole_G.number_of_edges()}")

# --- data exports (robust) ---
pld.assign_node_colors_for_gml(whole_G)
nx.write_gml(whole_G, f"{OUTDIR}/whole_network2.gml")
try:
    pld.write_filtered_gml_abundant_phylum(whole_G, f"{OUTDIR}/whole_network2_abundant_phylum.gml", min_host_per_phylum=5)
except Exception as e:
    print("  abundant-phylum gml skipped:", e)
pld.export_network_to_excel(whole_G, f"{OUTDIR}/network.xlsx")
pld.export_network_to_json(whole_G, f"{OUTDIR}/network.json")
try:
    pld.profile_network(whole_G, ctg_taxa_dict, clu_mge_dict)
except Exception as e:
    print("  profile_network skipped:", e)

gc_df = pd.DataFrame(gc_data, columns=["MGE", "host", "MGE_type", "MGE_gc", "host_gc",
        "cos_sim", "MGE_cov", "host_cov", "environment", "sample", "mge_len", "host_taxa", "host_phylum"])
gc_df = gc_df.sort_values(by="MGE_cov", ascending=False)
gc_df.to_csv(f"{OUTDIR}/mge_host_gc_cov.csv", index=False)
print(f"  wrote mge_host_gc_cov.csv ({len(gc_df)} linkages)")

# --- component analysis (largest component + binary matrix) ---
try:
    ccs = sorted(nx.weakly_connected_components(whole_G), key=len, reverse=True)
    print(f"  weakly-connected components: {len(ccs)}; largest={len(ccs[0]) if ccs else 0}")
    pld.component_analysis(ccs, all_host_clu_lineage_dict, ctg_taxa_dict, OUTDIR, clu_host_dict)
except Exception as e:
    print("  component_analysis skipped:", e)

# --- optional plots (skip on missing graphviz/seaborn deps) ---
for fn, args in [(pld.plot_network2, (whole_G, OUTDIR)), (pld.plot_gc, (gc_df, OUTDIR)),
                 (pld.analyze_MGEs, (gc_df, mge_clu_dict))]:
    try:
        fn(*args)
    except Exception as e:
        print(f"  {fn.__name__} skipped:", e)

print("DONE network rebuild ->", OUTDIR)
