import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite

def get_edge(host_sum_file, MGE_type_dict, bin2anno_dict):
    """
    Get the edge data from the host summary file.
    """
    # Read the host summary file
    host_sum = pd.read_csv(host_sum_file)
    
    G = nx.Graph()
    ## get node file
    for index, row in host_sum.iterrows():
        if row["pvalue"] >= 0.05:
            continue
        if row['MGE'] in MGE_type_dict:
            G.add_node(row['MGE'], label=row['MGE'], type=MGE_type_dict[row['MGE']])
        if row['host'] in bin2anno_dict:
            G.add_node(row['host'], label=row['host'], type=bin2anno_dict[row['host']])

        G.add_edge(row['MGE'], row['host'], weight=row['final_score'], type='share')
    plot(G)

def plot(G):
    import matplotlib.patches as mpatches

    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    # Use a large canvas
    plt.figure(figsize=(16, 12))  # Increased canvas size

    # pos = nx.spring_layout(G, k=0.15, iterations=20)
    # pos = nx.kamada_kawai_layout(G)


    # Separate nodes by type
    virus_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'virus']
    plasmid_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'plasmid']
    novel_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'novel']
    host_nodes = [n for n, d in G.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel']]


    # ## use bipartite layout, where hosts are on one side and MGEs on the other
    # pos = {}
    # # Assign positions for virus and plasmid nodes
    # for i, node in enumerate(virus_nodes + plasmid_nodes):
    #     pos[node] = (0, i * 0.1)
    # # Assign positions for host nodes
    # for i, node in enumerate(host_nodes):
    #     pos[node] = (2, i * 0.1)

    # Assume virus_nodes + plasmid_nodes are one set, host_nodes are the other
    top_nodes = virus_nodes + plasmid_nodes + novel_nodes
    pos = nx.bipartite_layout(G, top_nodes)  

    # Assign a color to each host type (phylum)
    host_types = list(set([G.nodes[n]['type'] for n in host_nodes]))
    color_map = plt.get_cmap('tab20')
    host_type_to_color = {t: color_map(i % 20) for i, t in enumerate(host_types)}
    host_colors = [host_type_to_color[G.nodes[n]['type']] for n in host_nodes]

    # Draw edges
    # nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.5)
    edge_type_to_color = {
        'share': 'green',
        'unique': 'gray',
        'transduction': 'purple',
        # add more types and colors as needed
        None: 'gray'  # default
    }
    # Group edges by type
    edge_type_dict = {}
    for u, v, d in G.edges(data=True):
        etype = d.get('type', None)
        edge_type_dict.setdefault(etype, []).append((u, v))

    for etype, edges in edge_type_dict.items():
        nx.draw_networkx_edges(
            G, pos,
            edgelist=edges,
            edge_color=edge_type_to_color.get(etype, 'gray'),
            alpha=0.7,
            width=2,
            label=etype
        )


    # Draw nodes: viruses and plasmids as squares, hosts as circles colored by type
    nx.draw_networkx_nodes(G, pos, nodelist=virus_nodes, node_shape='s', node_color='red', label='Virus', node_size=100)
    nx.draw_networkx_nodes(G, pos, nodelist=plasmid_nodes, node_shape='s', node_color='blue', label='Plasmid', node_size=100)
    nx.draw_networkx_nodes(G, pos, nodelist=novel_nodes, node_shape='s', node_color='orange', label='Novel', node_size=100)
    nx.draw_networkx_nodes(G, pos, nodelist=host_nodes, node_shape='o', node_color=host_colors, label='Host', node_size=100)

    # Do not show labels

    # Legend for host types
    handles = [
        mpatches.Patch(color='red', label='Virus'),
        mpatches.Patch(color='blue', label='Plasmid')
    ] + [
        mpatches.Patch(color=host_type_to_color[t], label=t) for t in host_types
    ]
    plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.title("Network of Host-MGE Interactions")
    plt.tight_layout(pad=3.0)
    plt.savefig("../../tmp/results/network_test_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
def get_MGE_type(mge_file):
    """
    Get the MGE type from the plasmid and virus summary files.
    """
    plasmid_types = {}
    plasmid_sum = pd.read_csv(mge_file, sep='\t')

    for index, row in plasmid_sum.iterrows():
        plasmid_types[row['seq_name']] = row["type"]

    return plasmid_types

def read_gtdb(gatk):
    """
    Read the GTDB summary file and return a dictionary of contig to bin mapping.
    """
    gtdb_df = pd.read_csv(gatk, sep='\t')
    bin2anno_dict = {}
    for index, row in gtdb_df.iterrows():
        if row['user_genome'] != 'user_genome':
            anno = row['classification']
            # print (f"anno: {anno}")
            if re.search('Unclassified', anno):
                phylum = "Unclassified"
            else:
                phylum = anno.split(';')[1].strip()
            bin2anno_dict[row['user_genome']] = phylum
    return bin2anno_dict

if __name__ == "__main__":  
    prefix = "infant_1"
    work_dir = f"/home/shuaiw/borg/paper/run2/{prefix}/"
    mge_file = os.path.join(work_dir, "all_mge.tsv")
    gtdk_bac_file = os.path.join(work_dir, "GTDB/gtdbtk.bac120.summary.tsv")
    gtdk_arc_file = os.path.join(work_dir, "GTDB/gtdbtk.ar122.summary.tsv")
    gtdk_all_file = os.path.join(work_dir, "GTDB/gtdbtk.all.summary.tsv")
    host_summary_file = os.path.join(work_dir, f"{prefix}_methylation/host_summary.csv")
    if os.path.exists(gtdk_arc_file):
        ## merge gtdk_bac_file and gtdk_arc_file to gtdk_all_file
        os.system(
            f"cat {gtdk_bac_file} {gtdk_arc_file} > {gtdk_all_file}"
        )
    else:
        os.system(
            f"cp {gtdk_bac_file} {gtdk_all_file}"
        )
    bin2anno_dict = read_gtdb(gtdk_all_file)
    MGE_type_dict = get_MGE_type(mge_file)
    get_edge(host_summary_file, MGE_type_dict, bin2anno_dict)