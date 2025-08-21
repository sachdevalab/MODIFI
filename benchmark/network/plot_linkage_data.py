import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
import plotly.graph_objects as go


def get_edge(host_sum_file, MGE_type_dict, bin2anno_dict, gc_data, environment, prefix):
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
        if row['final_score'] <= 0.6:
            continue
        gc_data.append([row['MGE_gc'], row['host_gc'], row['cos_sim'], row['MGE_cov'], row['host_cov'], environment, prefix])
        both_nodes = 0
        if row['host'] not in bin2anno_dict:
            continue
        if bin2anno_dict[row['host']] == "Unclassified":
            continue
        if row['MGE'] in MGE_type_dict:
            G.add_node(row['MGE'], label=row['MGE'], type=MGE_type_dict[row['MGE']])
            both_nodes += 1
        if row['host'] in bin2anno_dict:
            G.add_node(row['host'], label=row['host'], type=bin2anno_dict[row['host']])
            both_nodes += 1
        if both_nodes == 2:
            G.add_edge(row['MGE'], row['host'], weight=row['final_score'], type='share')

    return G, gc_data
    
def plot_network(G):
    import matplotlib.patches as mpatches
    
    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    
    # Use a large canvas for better visualization
    plt.figure(figsize=(16, 12))
    
    # Separate nodes by type
    virus_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'virus']
    plasmid_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'plasmid']
    novel_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'novel']
    host_nodes = [n for n, d in G.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel']]
    
    # Use spring layout with more iterations and fixed seed for reproducibility
    pos = nx.spring_layout(G, k=0.3, iterations=200, seed=42)
    
    # Assign a color to each host type (phylum)
    host_types = list(set([G.nodes[n]['type'] for n in host_nodes]))
    color_map = plt.get_cmap('tab20')
    host_type_to_color = {t: color_map(i % 20) for i, t in enumerate(host_types)}
    host_colors = [host_type_to_color[G.nodes[n]['type']] for n in host_nodes]

    # Draw edges first (so they appear behind nodes)
    nx.draw_networkx_edges(G, pos, edge_color='black', alpha=0.3, width=0.5)
    
    # Draw nodes with larger sizes and better visibility
    node_size = 50  # Larger nodes
    if virus_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=virus_nodes, node_shape='s', 
                              node_color='red', node_size=node_size, alpha=0.8)
    if plasmid_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=plasmid_nodes, node_shape='s', 
                              node_color='blue', node_size=node_size, alpha=0.8)
    if novel_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=novel_nodes, node_shape='s', 
                              node_color='black', node_size=node_size, alpha=0.8)
    if host_nodes:
        nx.draw_networkx_nodes(G, pos, nodelist=host_nodes, node_shape='o', 
                              node_color=host_colors, node_size=node_size, alpha=0.8)
    
    # Create legend
    legend_elements = []
    if virus_nodes:
        legend_elements.append(mpatches.Patch(color='red', label=f'Virus ({len(virus_nodes)})'))
    if plasmid_nodes:
        legend_elements.append(mpatches.Patch(color='blue', label=f'Plasmid ({len(plasmid_nodes)})'))
    if novel_nodes:
        legend_elements.append(mpatches.Patch(color='black', label=f'Novel ({len(novel_nodes)})'))
    
    # Add host type legend (limit to top 10 to avoid clutter)
    sorted_host_types = sorted(host_types, key=lambda t: sum(1 for n in host_nodes if G.nodes[n]['type'] == t), reverse=True)
    for i, t in enumerate(sorted_host_types[:10]):  # Show only top 10 host types
        count = sum(1 for n in host_nodes if G.nodes[n]['type'] == t)
        legend_elements.append(mpatches.Patch(color=host_type_to_color[t], label=f'{t} ({count})'))
    
    if len(sorted_host_types) > 10:
        legend_elements.append(mpatches.Patch(color='lightgray', label=f'... +{len(sorted_host_types)-10} more'))
    
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', 
              borderaxespad=0., fontsize=10)
    
    plt.axis('off')
    plt.title("Network of Host-MGE Interactions", fontsize=16, pad=20)
    plt.tight_layout()
    plt.savefig("../../tmp/results/network_test_plot.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot(G):
    import matplotlib.patches as mpatches

    print(f"Number of nodes: {G.number_of_nodes()}")
    print(f"Number of edges: {G.number_of_edges()}")
    # Use a large canvas
    plt.figure(figsize=(10, 8))  # Increased canvas size


    # Separate nodes by type
    virus_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'virus']
    plasmid_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'plasmid']
    novel_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'novel']
    host_nodes = [n for n, d in G.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel']]



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
        'share': 'gray',
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
            width=1,
            label=etype
        )


    # Draw nodes: viruses and plasmids as squares, hosts as circles colored by type
    nx.draw_networkx_nodes(G, pos, nodelist=virus_nodes, node_shape='s', node_color='red', label='Virus', node_size=100)
    nx.draw_networkx_nodes(G, pos, nodelist=plasmid_nodes, node_shape='s', node_color='blue', label='Plasmid', node_size=100)
    nx.draw_networkx_nodes(G, pos, nodelist=novel_nodes, node_shape='s', node_color='black', label='Novel', node_size=100)
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

def plot_gc(df):
    ## sort df by sample
    df = df.sort_values(by='sample')
    ## using seaborn to plot three subplots, first is MGE_gc vs host_gc, second is MGE_cov vs host_cov, third is a box plot show cos_sim in each environment
    ## use set2 as color palette
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set_palette("Set2")
    fig, axs = plt.subplots(3, 1, figsize=(6, 19))
    sns.scatterplot(data=df, x='MGE_gc', y='host_gc', ax=axs[0], hue='environment')
    # Set same limits for x and y axes
    min_gc = min(df['MGE_gc'].min(), df['host_gc'].min())
    max_gc = max(df['MGE_gc'].max(), df['host_gc'].max())
    axs[0].set_xlim(min_gc, max_gc)
    axs[0].set_ylim(min_gc, max_gc)
    axs[0].set_xlabel('MGE GC content')
    axs[0].set_ylabel('Host GC content')
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3, fontsize=9)

    sns.scatterplot(data=df, x='MGE_cov', y='host_cov', ax=axs[1], hue='environment')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    # axs[1].set_title('MGE Coverage vs Host Coverage')
    axs[1].set_xlabel('MGE Coverage (log scale)')
    axs[1].set_ylabel('Host Coverage (log scale)')
    axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3, fontsize=9)

    sns.boxplot(data=df, x='sample', y='cos_sim', ax=axs[2], hue='environment')
    # axs[2].set_title('Cosine Similarity Distribution by Environment')
    axs[2].set_xlabel('Environment')
    axs[2].set_ylabel('Cosine Similarity')
    axs[2].tick_params(axis='x', rotation=90)
    axs[2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=3, fontsize=9)

    plt.tight_layout()
    plt.savefig('../../tmp/results/mge_host_gc_content.pdf')
    plt.close()

if __name__ == "__main__":  
    meta_file = "/home/shuaiw/Methy/assembly_pipe/prefix_table.tab"
    sample_env_dict = read_metadata(meta_file)

    prefix = "96plex"
    ## the merged graph 
    whole_G = nx.Graph()
    gc_data = []
    # for prefix in ["96plex", "infant_1", "SRR23446539_sugarcane","ERR12723528_mice"]:
    all_dir = "/home/shuaiw/borg/paper/run2/"
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        # print (f"Processing {prefix}...")
        # work_dir = f"{all_dir}/{prefix}/{prefix}_methylation2"

        print(f"Processing {prefix}...")
        work_dir = f"/home/shuaiw/borg/paper/run2/{prefix}/"
        mge_file = os.path.join(work_dir, "all_mge.tsv")
        ## check if the mge_file exists
        if not os.path.exists(mge_file):
            print(f"File {mge_file} does not exist. Skipping {prefix}.")
            continue
        gtdk_bac_file = os.path.join(work_dir, "GTDB/gtdbtk.bac120.summary.tsv")
        gtdk_arc_file = os.path.join(work_dir, "GTDB/gtdbtk.ar122.summary.tsv")
        gtdk_all_file = os.path.join(work_dir, "GTDB/gtdbtk.all.summary.tsv")
        host_summary_file = os.path.join(work_dir, f"{prefix}_methylation2/host_summary.csv")
        print (f"Host summary file: {host_summary_file}")
        ## skip if the host_summary_file is empty, count the number of lines
        ## host sum lines
        line_num = sum(1 for line in open(host_summary_file) if line.strip())
        if line_num < 2:
            print(f"Host summary file {host_summary_file} is empty or has only header. Skipping {prefix}.")
            continue

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
        # print (bin2anno_dict)
        MGE_type_dict = get_MGE_type(mge_file)
        G, gc_data = get_edge(host_summary_file, MGE_type_dict, bin2anno_dict, gc_data, sample_env_dict[prefix], prefix)
        whole_G = nx.compose(whole_G, G)
        # break
    ## print the node with top 10 degree
    top_nodes = sorted(whole_G.degree, key=lambda x: x[1], reverse=True)[:10]
    print("Top 10 nodes by degree:")
    for node, degree in top_nodes:
        print(f"{node}: {degree}")
    plot_network(whole_G)

    gc_df = pd.DataFrame(gc_data, columns=["MGE_gc", "host_gc", "cos_sim", "MGE_cov", "host_cov", "environment", "sample"])
    plot_gc(gc_df)