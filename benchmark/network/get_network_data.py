import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt



def get_edge(plasmid_types={}, virus_types={}, bin2anno_dict={}, cutoff=0.45):
    """
    Get the edge data from the host summary file.
    """
    # Read the host summary file
    host_sum = pd.read_csv(host_sum_file)
    
    G = nx.Graph()
    ## get node file
    for index, row in host_sum.iterrows():
        if row['final_score'] < cutoff:
            continue
        if row['MGE'] in plasmid_types:
            G.add_node(row['MGE'], label=row['MGE'], type="Plasmid")
        elif row['MGE'] in virus_types:
            G.add_node(row['MGE'], label=row['MGE'], type="Virus")
        else:
            G.add_node(row['MGE'], label=row['MGE'], type="Virus")
        if row['host'] in bin2anno_dict:

            G.add_node(row['host'], label=row['host'], type=bin2anno_dict[row['host']])
        else:
            G.add_node(row['host'], label=row['host'], type="Unknown")

        if row['both_link'] == 1:
            G.add_edge(row['MGE'], row['host'], weight=row['final_score'], type='share')
        elif row['both_link'] == 0:
            G.add_edge(row['MGE'], row['host'], weight=row['final_score'], type='unique')
    plot(G)
    nx.write_gexf(G,gexf)

def plot(G):
    import matplotlib.patches as mpatches

    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    # Use a large canvas
    plt.figure(figsize=(16, 12))  # Increased canvas size

    # pos = nx.spring_layout(G, k=0.15, iterations=20)
    pos = nx.kamada_kawai_layout(G)


    # Separate nodes by type
    virus_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'Virus']
    plasmid_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'Plasmid']
    host_nodes = [n for n, d in G.nodes(data=True) if d.get('type') not in ['Virus', 'Plasmid', 'MGE']]

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
    plt.savefig("../../tmp/results/network_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    

def get_MGE_type(plasmid_sum, virus_sum):
    """
    Get the MGE type from the plasmid and virus summary files.
    """
    plasmid_types = {}
    virus_types = {}
    plasmid_sum = pd.read_csv(plasmid_sum, sep='\t')
    virus_sum = pd.read_csv(virus_sum, sep='\t')

    for index, row in plasmid_sum.iterrows():
        plasmid_types[row['seq_name']] = 'plasmid'
    for index, row in virus_sum.iterrows():
        virus_types[row['seq_name']] = 'virus'

    return plasmid_types, virus_types

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

def classify_cow_MGE():
    MGE_type_file = "/home/shuaiw/borg/pengfan/contigs/All_Circular_Elements_genomad_summary.tsv"
    plasmid_types={}
    virus_types={}
    for line in open(MGE_type_file, "r"):
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        if len(line) < 3:
            continue
        MGE_name = line[0]
        MGE_type = line[-1]
        if MGE_type == "plasmid":
            plasmid_types[MGE_name] = MGE_type
        elif MGE_type == "virus":
            virus_types[MGE_name] = MGE_type
    return plasmid_types, virus_types

if __name__ == "__main__":  
    # host_sum_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/host_summary.csv"
    # gexf = "/home/shuaiw/borg/paper/network/RuReacBro_20230708_11_72h_20_bin.gexf"

    # host_sum_file = "/home/shuaiw/borg/bench/soil/run1/host_summary.csv"
    # gexf = "/home/shuaiw/borg/paper/network/soil_run1_bin.gexf"
    # plasmid_sum = "/home/shuaiw/methylation/data/borg/contigs/genomad/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_summary/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_plasmid_summary.tsv"
    # virus_sum = "/home/shuaiw/methylation/data/borg/contigs/genomad/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_summary/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_virus_summary.tsv"
    # gatk = "/home/shuaiw/methylation/data/borg/contigs/GTDB/gtdbtk.all.summary.tsv"
    # plasmid_types, virus_types = get_MGE_type(plasmid_sum, virus_sum)
    # bin2anno_dict = read_gtdb(gatk)
    # get_edge(plasmid_types, virus_types, bin2anno_dict)  

    # host_sum_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/host_summary.csv"
    # host_sum_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_12_72h_200ppm_r2_HMW_LR_bin/host_summary.csv"
    host_sum_file = "/home/shuaiw/methylation/data/borg/pengfan/total_summary_compare.csv"
    gexf = "/home/shuaiw/borg/paper/network/RuReacBro_20230708_11_72h_20_bin.gexf"
    gatk = "/home/shuaiw/borg/pengfan/contigs/gatk_all.summary.tsv"
    bin2anno_dict = read_gtdb(gatk)
    plasmid_types, virus_types = classify_cow_MGE()
    get_edge(plasmid_types, virus_types, bin2anno_dict=bin2anno_dict, cutoff=0.45)  