import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import plotly.graph_objects as go
from collections import defaultdict
import matplotlib.patches as mpatches
import random
import sys
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_detail_taxa_name,get_unique_motifs, My_sample, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa



def get_edge(cluster_anno_dict, MGE_type_dict, gc_data, environment, prefix, 
             host_clu_dict, mge_clu_dict, sample_obj, ctg_taxa_dict):
    """
    Get the edge data from the host summary file.
    """
    our_linkages, our_ctg_linkages, linkage_info_list = sample_obj.read_linkage_dict()

    
    G = nx.Graph()
    for linkage_obj in linkage_info_list:

        host_lineage = ctg_taxa_dict[linkage_obj.host] if linkage_obj.host in ctg_taxa_dict else "NA"
        host_taxa = get_detail_taxa_name(host_lineage)

        gc_data.append([linkage_obj.mge, linkage_obj.host, linkage_obj.MGE_gc, linkage_obj.host_gc, linkage_obj.cos_sim, 
                        linkage_obj.MGE_cov, linkage_obj.host_cov, environment, prefix, linkage_obj.mge_len, host_taxa])

        if linkage_obj.host in host_clu_dict:
            host_clu = host_clu_dict[linkage_obj.host]
        else:
            host_clu = linkage_obj.host
        if linkage_obj.mge in MGE_type_dict:
            # G.add_node(row['MGE'], label=row['MGE'], type=MGE_type_dict[row['MGE']])
            G.add_node(mge_clu_dict[linkage_obj.mge], label=linkage_obj.mge, type=MGE_type_dict[linkage_obj.mge])

        represent_ctg_lineage = ctg_taxa_dict[linkage_obj.host] if linkage_obj.host in ctg_taxa_dict else "NA"
        ctg_phylum = classify_taxa(represent_ctg_lineage, "phylum")
        G.add_node(host_clu, label=linkage_obj.host, type=ctg_phylum)
        cluster_anno_dict[host_clu] = ctg_phylum
        G.add_edge(mge_clu_dict[linkage_obj.mge], host_clu, weight=linkage_obj.final_score, type='share')

    return G, gc_data, cluster_anno_dict
    
def plot_network2(G,paper_fig_dir ):

    
    plt.figure(1, figsize=(10, 12))
    
    # layout graphs with positions using graphviz neato
    pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
    
    # Separate nodes by type
    virus_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'virus']
    plasmid_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'plasmid']
    novel_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'novel']
    host_nodes = [n for n, d in G.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel']]
    
    # Draw edges first (so they appear behind nodes)
    nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.8, width=0.5)
    
    # Define colors for each node type
    node_size = 50
    
    # Assign a unique color to each individual host cluster
    unique_host_clusters = list(set(host_nodes))  # Get unique host cluster IDs
    color_map = plt.get_cmap('tab20')
    host_cluster_to_color = {cluster: color_map(i % 20) for i, cluster in enumerate(unique_host_clusters)}
    host_colors = [host_cluster_to_color[n] for n in host_nodes]
    
    # Group clusters by phylum for legend
    host_types = list(set([G.nodes[n]['type'] for n in host_nodes]))
    host_type_to_clusters = {t: [n for n in host_nodes if G.nodes[n]['type'] == t] for t in host_types}
    
    # Draw different node types with different shapes and colors
    if host_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=host_nodes, node_shape='s',  # square
                              node_color=host_colors, node_size=node_size, alpha=0.8, label='Host')
    
    if plasmid_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=plasmid_nodes, node_shape='o',  # circle
                              node_color='blue', node_size=node_size, alpha=0.8, label='Plasmid')
    
    if virus_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=virus_nodes, node_shape='h',  # hexagon
                              node_color='red', node_size=node_size, alpha=0.8, label='Virus')
    
    if novel_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=novel_nodes, node_shape='^',  # triangle
                              node_color='orange', node_size=node_size, alpha=0.8, label='Novel')
    
    # Create legend with shape indicators
    from matplotlib.lines import Line2D
    legend_elements = []
    
    if plasmid_nodes:
        legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', 
                                    markersize=8, label=f'Plasmid ({len(plasmid_nodes)})'))
    if virus_nodes:
        legend_elements.append(Line2D([0], [0], marker='h', color='w', markerfacecolor='red', 
                                    markersize=8, label=f'Virus ({len(virus_nodes)})'))
    if novel_nodes:
        legend_elements.append(Line2D([0], [0], marker='^', color='w', markerfacecolor='orange', 
                                    markersize=8, label=f'Novel ({len(novel_nodes)})'))
    
    # Add host type legend (limit to top 10 to avoid clutter)
    if host_nodes:
        sorted_host_types = sorted(host_types, key=lambda t: sum(1 for n in host_nodes if G.nodes[n]['type'] == t), reverse=True)
        for i, t in enumerate(sorted_host_types[:6]):  # Show only top 10 host types
            count = sum(1 for n in host_nodes if G.nodes[n]['type'] == t)
            # Use the color of the first cluster of this phylum type as representative color
            representative_cluster = host_type_to_clusters[t][0]
            legend_elements.append(Line2D([0], [0], marker='s', color='w', markerfacecolor=host_cluster_to_color[representative_cluster], 
                                        markersize=8, label=f'{t} ({count})'))
        
        if len(sorted_host_types) > 6:
            legend_elements.append(Line2D([0], [0], marker='s', color='w', markerfacecolor='lightgray', 
                                        markersize=8, label=f'... +{len(sorted_host_types)-10} more host types'))
    
    if legend_elements:
        plt.legend(handles=legend_elements, bbox_to_anchor=(0.5, -0.05), loc='upper center', 
                  borderaxespad=0., fontsize=10, ncol=3)
    
    plt.title(f"Host-MGE Network\nNodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}", 
              fontsize=14, pad=20)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(f"{paper_fig_dir}/network_test_plot2.png", dpi=300, bbox_inches='tight')
    plt.close()

def read_metadata(meta_file):
    sample_env_dict = {}
    with open(meta_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            sample = parts[1]
            env = parts[3]
            sample_env_dict[sample] = env
    return sample_env_dict

def plot_gc(df, paper_fig_dir):
    ## sort df by sample
    df = df.sort_values(by='sample')
    ## using seaborn to plot three subplots, first is MGE_gc vs host_gc, second is MGE_cov vs host_cov, third is a box plot show cos_sim in each environment
    ## use set2 as color palette

    sns.set_palette("Set2")
    fig, axs = plt.subplots(1, 3, figsize=(16, 6))
    sns.scatterplot(data=df, x='MGE_gc', y='host_gc', ax=axs[0], hue='environment')
    # Set same limits for x and y axes
    min_gc = min(df['MGE_gc'].min(), df['host_gc'].min())
    max_gc = max(df['MGE_gc'].max(), df['host_gc'].max())
    axs[0].set_xlim(min_gc, max_gc)
    axs[0].set_ylim(min_gc, max_gc)
    axs[0].set_xlabel('MGE GC content')
    axs[0].set_ylabel('Host GC content')
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, fontsize=9)

    sns.scatterplot(data=df, x='MGE_cov', y='host_cov', ax=axs[1], hue='environment', legend=False)
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    # axs[1].set_title('MGE Coverage vs Host Coverage')
    axs[1].set_xlabel('MGE Coverage (log scale)')
    axs[1].set_ylabel('Host Coverage (log scale)')


    # Calculate median cos_sim for each environment and sort
    env_order = df.groupby('environment')['cos_sim'].median().sort_values(ascending=False).index
    sns.boxplot(data=df, x='environment', y='cos_sim', ax=axs[2], hue='environment', legend=False, order=env_order)
    # axs[2].set_title('Cosine Similarity Distribution by Environment')
    axs[2].set_xlabel('Environment')
    axs[2].set_ylabel('Cosine Similarity')
    axs[2].tick_params(axis='x', rotation=45, labelsize=9)

    plt.tight_layout()
    plt.savefig(f"{paper_fig_dir}/mge_host_gc_content.png", dpi=300)
    plt.close()

def read_drep_cluster(drep_clu_file):
    
    drep_clu_dict = defaultdict(list)
    host_clu_dict = {}
    df = pd.read_csv(drep_clu_file)
    for index, row in df.iterrows():
        # drep_clu_dict[row['genome'][:-3]] = row['secondary_cluster']
        contig = row['genome'][:-3]
        drep_clu_dict[row['secondary_cluster']].append(row['genome'][:-3])
        host_clu_dict[contig] = row['secondary_cluster']
    return drep_clu_dict, host_clu_dict

def read_mge_cluster(mge_clu_file):
    mge_clu_dict = {}
    df = pd.read_csv(mge_clu_file, sep="\t", header=None)
    for index, row in df.iterrows():
        cluster = row[0]
        mges = row[1].split(",")
        for mge in mges:
            mge_clu_dict[mge] = cluster
    return mge_clu_dict

def count_cross_phylum(whole_G):
    ## analyze the node in ['virus', 'plasmid', 'novel'], count how many of them are linked to multiple phyla for each of them
    ## also count the number of nodes in each type in  these three types
    ## also count the number of nodes number of nodes in each type in these three types
    cross_phylum_count = {'virus': 0, 'plasmid': 0, 'novel': 0}
    type_count = {'virus': 0, 'plasmid': 0, 'novel': 0}
    for node, degree in whole_G.degree:
        if whole_G.nodes[node]['type'] in ['virus', 'plasmid', 'novel']:
            type_count[whole_G.nodes[node]['type']] += 1
            neighbors = list(whole_G.neighbors(node))
            phyla = set()
            for neighbor in neighbors:
                if whole_G.nodes[neighbor]['type'] not in ['virus', 'plasmid', 'novel']:
                    phyla.add(whole_G.nodes[neighbor]['type'])
            if len(phyla) > 1:
                print (f"{node} ({whole_G.nodes[node]['type']}) is linked to multiple phyla: {phyla}")
                cross_phylum_count[whole_G.nodes[node]['type']] += 1
    print("Cross-phylum linked MGEs:")

    print(cross_phylum_count)
    print("Number of nodes in each type:")
    print(type_count)

def profile_network(whole_G, ctg_taxa_dict):
    ## print the node with top 10 degree
    top_nodes = sorted(whole_G.degree, key=lambda x: x[1], reverse=True)[:10]
    print("Top 10 nodes by degree:")
    for node, degree in top_nodes:
        node_type = whole_G.nodes[node].get('type', 'Unknown')
        if node_type in ['virus', 'plasmid', 'novel']:
            print(f"{node}: {degree} (MGE type: {node_type})")
        else:
            ## get the label of the node
            node_label = whole_G.nodes[node].get('label', 'Unknown')
            represent_ctg_lineage = ctg_taxa_dict[node_label] if node_label in ctg_taxa_dict else "NA"
            ctg_species = classify_taxa(represent_ctg_lineage, "species")
            ctg_taxa = get_detail_taxa_name(represent_ctg_lineage)
            print(f"{node} {node_label}: {degree} (Host annotation: {ctg_species}, {ctg_taxa})")
    
    ## see which host has linked the most virus
    # host_virus_link_count = defaultdict(int)
    # for node, degree in whole_G.degree:
    #     if whole_G.nodes[node]['type'] not in ['virus', 'plasmid', 'novel']:
    #         neighbors = list(whole_G.neighbors(node))
    #         for neighbor in neighbors:
    #             if whole_G.nodes[neighbor]['type'] == 'virus':
    #                 host_virus_link_count[node] += 1
    # # get the top 10 hosts that linked the most virus
    # top_hosts = sorted(host_virus_link_count.items(), key=lambda x: x[1], reverse=True)[:10]
    # print("Top 10 hosts linking the most viruses:")
    # for host, count in top_hosts:
    #     node_label = whole_G.nodes[host].get('label', 'Unknown')
    #     represent_ctg_lineage = ctg_taxa_dict[node_label] if node_label in ctg_taxa_dict else "NA"
    #     ctg_species = classify_taxa(represent_ctg_lineage, "species")
    #     ctg_genus = classify_taxa(represent_ctg_lineage, "genus")
    #     print(f"{host} {node_label}: linked to {count} viruses (Host annotation: {ctg_species}, {ctg_genus})")
    # count_cross_phylum(whole_G)
    # # ## print the MGE with degree > 1
    # # print("MGEs with degree > 1:")
    # for node, degree in whole_G.degree:
    #     if degree > 1 and whole_G.nodes[node]['type'] in ['virus']:
    #         print(f"virus {node} ({whole_G.nodes[node]['type']}): {degree}")
    #     # if degree > 1 and whole_G.nodes[node]['type'] in ['plasmid']:
    #     #     print(f"plasmid {node} ({whole_G.nodes[node]['type']}): {degree}")
    #     # if degree > 1 and whole_G.nodes[node]['type'] in ['novel']:
    #     #     print(f"novel {node} ({whole_G.nodes[node]['type']}): {degree}")
    
    #         ### print the linked nodes to infant_15_839_L
    #         target_node = node #"cow_bioreactor_2_2972_L"
    #         if target_node in whole_G:
    #             neighbors = list(whole_G.neighbors(target_node))
    #             print(f"Neighbors of {target_node}: {neighbors}")
    #         print ("#########################")

if __name__ == "__main__":  
    ANI = 99
    meta_file = "/home/shuaiw/mGlu/assembly_pipe/prefix_table.tab"
    mge_clu_file = "/home/shuaiw/borg/paper/MGE/cluster/megablast.cluster.95ani.tsv"
    paper_fig_dir = f"../../tmp/figures/multi_env_linkage/network_{ANI}/"
    drep_clu_file = f"/home/shuaiw/borg/paper/specificity/dRep_{ANI}_out/data_tables/Cdb.csv"
    all_dir = "/home/shuaiw/borg/paper/run2/"

    ## mkdir paper_fig_dir if not exists
    if not os.path.exists(paper_fig_dir):
        os.makedirs(paper_fig_dir)

    sample_env_dict = read_metadata(meta_file)
    drep_clu_dict, host_clu_dict = read_drep_cluster(drep_clu_file)
    mge_clu_dict = read_mge_cluster(mge_clu_file)
    ctg_taxa_dict = get_ctg_taxa(all_dir)


    ## the merged graph 
    whole_G = nx.Graph()
    gc_data = []
    cluster_anno_dict = {}
    
    for my_dir in os.listdir(all_dir):
        prefix = my_dir

        if prefix in [ "ERR5621427_sludge", "ERR5621429_sludge", "ERR5621430_sludge"]:
            continue

        sample_obj = My_sample(prefix, all_dir)
        print(f"Processing {prefix}...")


        MGE_type_dict = sample_obj.read_MGE()
        G, gc_data, cluster_anno_dict = get_edge(cluster_anno_dict, MGE_type_dict, gc_data, 
                                                 sample_env_dict[prefix], prefix, host_clu_dict, 
                                                 mge_clu_dict, sample_obj, ctg_taxa_dict)
        whole_G = nx.compose(whole_G, G)
        # break



    plot_network2(whole_G, paper_fig_dir)
    profile_network(whole_G, ctg_taxa_dict)

    gc_df = pd.DataFrame(gc_data, columns=["MGE", "host", "MGE_gc", "host_gc", 
                                           "cos_sim", "MGE_cov", "host_cov", 
                                           "environment", "sample", "mge_len", "host_taxa"])
    ## sort gc_df by mge_len descending
    gc_df = gc_df.sort_values(by='mge_len', ascending=False)
    ## save gc_df to a csv file
    gc_df.to_csv(f"{paper_fig_dir}/mge_host_gc_cov.csv", index=False)
    plot_gc(gc_df, paper_fig_dir)

    