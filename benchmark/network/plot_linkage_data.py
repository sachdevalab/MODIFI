import json
import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from networkx.algorithms import bipartite
import plotly.graph_objects as go
from collections import defaultdict
import matplotlib.patches as mpatches
import random
import sys
import seaborn as sns
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_detail_taxa_name,get_unique_motifs, My_sample, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'specificity'))
from profile_good_ctgs import BAC_PHYLUM_COLORS, ARC_PHYLUM_COLORS

def export_mge_degrees(G, out_csv):
    """
    Count the degree of each virus and plasmid node in the graph G and output to a CSV file.
    Only virus and plasmid nodes are included.
    """
    rows = []
    for node, degree in G.degree:
        node_type = G.nodes[node].get('type', None)
        if node_type in ['virus', 'plasmid']:
            rows.append({'node': node, 'type': node_type, 'degree': int(degree)})
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(out_csv, index=False)
    else:
        print('No virus or plasmid nodes found in the graph.')


def export_network_to_excel(G, out_path):
    """
    Export network to an Excel file with two sheets: Nodes and Edges.
    Nodes: Node_ID, Label, Type
    Edges: Source, Target, Weight, Type
    """
    nodes_rows = []
    for node, data in G.nodes(data=True):
        nodes_rows.append({
            'Node_ID': str(node),
            'Label': data.get('label', node),
            'Type': data.get('type', ''),
        })
    nodes_df = pd.DataFrame(nodes_rows)

    edges_rows = []
    for u, v, data in G.edges(data=True):
        edges_rows.append({
            'Source': str(u),
            'Target': str(v),
            'Weight': data.get('weight', None),
            'Type': data.get('type', ''),
        })
    edges_df = pd.DataFrame(edges_rows)

    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        nodes_df.to_excel(writer, sheet_name='Nodes', index=False)
        edges_df.to_excel(writer, sheet_name='Edges', index=False)
    print(f"Network exported to {out_path} (Nodes: {len(nodes_df)}, Edges: {len(edges_df)})")


def export_network_to_json(G, out_path):
    """
    Export network to JSON in the format:
    {"nodes": [{"id": ..., "attributes": {"Name": label, "Type": type}}, ...],
     "edges": [{"id": idx, "source": ..., "target": ..., "attributes": {"Weight": ..., "Type": ...}}, ...]}
    """
    nodes_list = []
    for node, data in G.nodes(data=True):
        label = data.get('label', node)
        node_type = data.get('type', '')
        nodes_list.append({
            "id": str(node),
            "attributes": {
                "Name": str(label),
                "Type": str(node_type),
            }
        })

    edges_list = []
    for idx, (u, v, data) in enumerate(G.edges(data=True)):
        attrs = {}
        if data.get('weight') is not None:
            attrs["Weight"] = float(data['weight'])
        if data.get('type') is not None:
            attrs["Type"] = str(data['type'])
        edges_list.append({
            "id": idx,
            "source": str(u),
            "target": str(v),
            "attributes": attrs,
        })

    out = {"nodes": nodes_list, "edges": edges_list}
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Network exported to {out_path} (Nodes: {len(nodes_list)}, Edges: {len(edges_list)})")


def _matplotlib_color_to_hex(color):
    """Convert matplotlib color (name or (r,g,b[,a]) 0-1) to hex string for GML/Gephi."""
    if isinstance(color, str) and color.startswith('#'):
        return color
    try:
        rgb = mcolors.to_rgb(color)
        return '#{:02X}{:02X}{:02X}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
    except Exception:
        return '#808080'


# MGE type colors for GML/Gephi (match legend: plasmid = salmon, virus = teal)
MGE_TYPE_COLORS = {
    'plasmid': "#F8766D",   # salmon pink / light coral
    'virus': "#00BFC4",     # teal / light sea green
}


def assign_node_colors_for_gml(G):
    """
    Set node color attributes on G so that GML export (and Gephi) can display them.
    Host nodes: phylum colors from profile_good_ctgs BAC_PHYLUM_COLORS and ARC_PHYLUM_COLORS.
    MGE nodes: plasmid/virus/novel colors (plasmid=salmon, virus=teal per legend).
    Sets 'graphics_fill' (hex) and 'color' (hex) on each node for Gephi compatibility.
    """
    virus_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'virus']
    plasmid_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'plasmid']
    host_nodes = [n for n, d in G.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel']]

    others_hex = BAC_PHYLUM_COLORS.get('Others', '#e6e6e6')
    host_phyla = set(G.nodes[n].get('type') for n in host_nodes)
    host_phylum_to_color = {}
    for p in host_phyla:
        p_norm = p.replace("p__", "") if p else p
        host_phylum_to_color[p] = (
            BAC_PHYLUM_COLORS.get(p_norm) or ARC_PHYLUM_COLORS.get(p_norm) or others_hex
        )

    for n in G.nodes:
        node_type = G.nodes[n].get('type')
        if n in virus_nodes:
            hexcolor = MGE_TYPE_COLORS['virus']
        elif n in plasmid_nodes:
            hexcolor = MGE_TYPE_COLORS['plasmid']

        elif node_type and node_type in host_phylum_to_color:
            hexcolor = host_phylum_to_color[node_type]
        else:
            hexcolor = others_hex
        G.nodes[n]['graphics_fill'] = hexcolor
        G.nodes[n]['color'] = hexcolor

    # Print only node categories present in this network (for manual Gephi Partition)
    type_to_color = {}
    for n in G.nodes:
        t = G.nodes[n].get('type')
        c = G.nodes[n].get('color')
        if t:
            type_to_color[t] = c
    mge_order = ['virus', 'plasmid', 'novel']
    keys_sorted = sorted(type_to_color.keys(),
                         key=lambda x: (0 if x in mge_order else 1, mge_order.index(x) if x in mge_order else 999, x))
    print("\n--- Node categories in network (Gephi Partition -> type; paste hex into color square) ---")
    for t in keys_sorted:
        print(f"  {t}\t{type_to_color[t]}")
    print("---\n")


def write_filtered_gml_abundant_phylum(G, out_path, min_host_per_phylum=5):
    """
    Write a filtered GML that keeps only host nodes in abundant phyla (> min_host_per_phylum nodes).
    MGE nodes are kept only if they have at least one edge to a kept host.
    """
    mge_types = {'virus', 'plasmid', 'novel'}
    # Host nodes: type = phylum (not MGE type)
    host_phylum_counts = defaultdict(int)
    for n, d in G.nodes(data=True):
        t = d.get('type')
        if t and t not in mge_types:
            host_phylum_counts[t] += 1
    abundant_phyla = {p for p, c in host_phylum_counts.items() if c > min_host_per_phylum}
    # Host nodes to keep: phylum in abundant_phyla
    host_nodes_keep = {n for n, d in G.nodes(data=True)
                       if d.get('type') not in mge_types and d.get('type') in abundant_phyla}
    # MGE nodes to keep: have at least one edge to a kept host (for DiGraph: successor in host_nodes_keep)
    mge_nodes_keep = set()
    for n, d in G.nodes(data=True):
        if d.get('type') not in mge_types:
            continue
        for succ in G.successors(n):
            if succ in host_nodes_keep:
                mge_nodes_keep.add(n)
                break
    nodes_to_keep = host_nodes_keep | mge_nodes_keep
    filtered_G = G.subgraph(nodes_to_keep).copy()
    nx.write_gml(filtered_G, out_path)
    print(f"Filtered GML (hosts only in phyla with >{min_host_per_phylum} nodes): "
          f"{filtered_G.number_of_nodes()} nodes, {filtered_G.number_of_edges()} edges -> {out_path}")


def collect_host_ctgs(prefix, sample_obj, ctg_taxa_dict,host_ctg_set):
    """
    Get the edge data from the host summary file.
    """
    # sample_obj.specificity_cutoff = 0.001
    # sample_obj.final_score_cutoff = 0.8
    our_linkages, our_ctg_linkages, linkage_info_list = sample_obj.read_linkage_dict()

    host_clu_lineage_dict = {}
    G = nx.Graph()
    for linkage_obj in linkage_info_list:

        host_lineage = ctg_taxa_dict[linkage_obj.host] if linkage_obj.host in ctg_taxa_dict else "NA"
        host_taxa = get_detail_taxa_name(host_lineage)
        ## get phylum of host
        host_phylum = classify_taxa(host_lineage, "phylum")
        host_lineage = ctg_taxa_dict[linkage_obj.host] if linkage_obj.host in ctg_taxa_dict else "NA"
        
        ## skip the linkage if host taxa is unclassfied
        if re.search("Unclassified", host_taxa):
            continue
        ctg_obj = My_contig(prefix, sample_obj.all_dir, linkage_obj.host)
        host_ctg_set.add(ctg_obj.ctg_ref)
    return host_ctg_set

def get_edge(cluster_anno_dict, MGE_type_dict, gc_data, environment, prefix, 
             host_clu_dict, mge_clu_dict, sample_obj, ctg_taxa_dict):
    """
    Get the edge data from the host summary file.
    Builds a directed graph with edges MGE -> host.
    """
    sample_obj.specificity_cutoff = 0.001
    sample_obj.final_score_cutoff = 0.8
    our_linkages, our_ctg_linkages, linkage_info_list = sample_obj.read_linkage_dict()

    host_clu_lineage_dict = {}
    G = nx.DiGraph()
    for linkage_obj in linkage_info_list:

        host_lineage = ctg_taxa_dict[linkage_obj.host] if linkage_obj.host in ctg_taxa_dict else "NA"
        host_taxa = get_detail_taxa_name(host_lineage)
        ## get phylum of host
        host_phylum = classify_taxa(host_lineage, "phylum")
        host_lineage = ctg_taxa_dict[linkage_obj.host] if linkage_obj.host in ctg_taxa_dict else "NA"
        
        ## skip the linkage if host taxa is unclassfied
        if re.search("Unclassified", host_taxa):
            continue
        # print (f"Processing linkage: host={linkage_obj.host}, mge={linkage_obj.mge} {MGE_type_dict[linkage_obj.mge]}")
        if MGE_type_dict[linkage_obj.mge] == "novel":
            continue
        if linkage_obj.host in host_clu_dict:
            host_clu = host_clu_dict[linkage_obj.host]
        else:
            host_clu = linkage_obj.host
        host_clu_lineage_dict[host_clu] = host_lineage
        if linkage_obj.mge in MGE_type_dict:
            # G.add_node(row['MGE'], label=row['MGE'], type=MGE_type_dict[row['MGE']])
            G.add_node(mge_clu_dict[linkage_obj.mge], label=linkage_obj.mge, type=MGE_type_dict[linkage_obj.mge])

        gc_data.append([linkage_obj.mge, linkage_obj.host, MGE_type_dict[linkage_obj.mge], linkage_obj.MGE_gc, linkage_obj.host_gc, linkage_obj.cos_sim, 
                        linkage_obj.MGE_cov, linkage_obj.host_cov, environment, prefix, linkage_obj.mge_len, host_taxa, host_phylum])

        represent_ctg_lineage = ctg_taxa_dict[linkage_obj.host] if linkage_obj.host in ctg_taxa_dict else "NA"
        ctg_phylum = classify_taxa(represent_ctg_lineage, "phylum")
        G.add_node(host_clu, label=linkage_obj.host, type=ctg_phylum)
        cluster_anno_dict[host_clu] = ctg_phylum
        # Directed edge: MGE -> host
        G.add_edge(mge_clu_dict[linkage_obj.mge], host_clu, weight=linkage_obj.final_score, type='share')

    return G, gc_data, cluster_anno_dict, host_clu_lineage_dict
    
def plot_network2(G, paper_fig_dir, fig_name="network_test_plot2"):

    
    # plt.figure(1, figsize=(12, 14))
    plt.figure(1, figsize=(16, 10))
    
    # layout graphs with positions using graphviz neato
    # pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
    pos = nx.nx_agraph.graphviz_layout(G, prog="fdp")
    
    # Separate nodes by type
    virus_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'virus']
    plasmid_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'plasmid']
    novel_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'novel']
    host_nodes = [n for n, d in G.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel']]
    
    # Draw edges first (so they appear behind nodes)
    nx.draw_networkx_edges(G, pos, edge_color='black', alpha=1, width=1)
    
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
    plt.savefig(f"{paper_fig_dir}/{fig_name}.png", dpi=300, bbox_inches='tight')
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
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    # use marker shape for MGE type and color for environment
    try:
        markers_map = {'plasmid': 'o', 'virus': 'X', 'novel': '^'}
        sns.scatterplot(data=df, x='MGE_gc', y='host_gc', ax=axs[0, 0], hue='environment', style='MGE_type', markers=markers_map)
    except Exception:
        sns.scatterplot(data=df, x='MGE_gc', y='host_gc', ax=axs[0, 0], hue='environment')
        # compute Pearson correlation for GC (MGE_gc vs host_gc) and annotate
    try:
        import numpy as _np
        from scipy.stats import pearsonr
        valid_gc = df[['MGE_gc', 'host_gc']].dropna()
        if len(valid_gc) > 1:
            r_gc, p_gc = pearsonr(valid_gc['MGE_gc'], valid_gc['host_gc'])
        else:
            r_gc, p_gc = _np.nan, _np.nan
    except Exception:
        try:
            import numpy as _np
            valid_gc = df[['MGE_gc', 'host_gc']].dropna()
            r_gc = _np.corrcoef(valid_gc['MGE_gc'].astype(float), valid_gc['host_gc'].astype(float))[0, 1] if len(valid_gc) > 1 else _np.nan
            p_gc = _np.nan
        except Exception:
            r_gc, p_gc = None, None
    try:
        if r_gc is not None:
            axs[0, 0].text(0.05, 0.95, f"r={r_gc:.2f}\n p={p_gc:.2g}", transform=axs[0, 0].transAxes, va='top', fontsize=9)
    except Exception:
        pass
    # regression lines removed per user request (keep only r/p annotation)
    # Set same limits for x and y axes
    min_gc = min(df['MGE_gc'].min(), df['host_gc'].min())
    max_gc = max(df['MGE_gc'].max(), df['host_gc'].max())
    axs[0, 0].set_xlim(min_gc, max_gc)
    axs[0, 0].set_ylim(min_gc, max_gc)
    axs[0, 0].set_xlabel('MGE GC content')
    axs[0, 0].set_ylabel('Host GC content')
    # move legend outside to the right to avoid overlap with points
    try:
        axs[0, 0].legend(loc='upper left', bbox_to_anchor=(1.02, 1), ncol=1, fontsize=9)
    except Exception:
        pass

    # use marker shape for MGE type and color for environment (keep legend off to avoid clutter)
    try:
        markers_map = {'plasmid': 'o', 'virus': 'X', 'novel': '^'}
        sns.scatterplot(data=df, x='MGE_cov', y='host_cov', ax=axs[0, 1], hue='environment', style='MGE_type', markers=markers_map, legend=False)
    except Exception:
        sns.scatterplot(data=df, x='MGE_cov', y='host_cov', ax=axs[0, 1], hue='environment', legend=False)
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_yscale('log')
    # axs[0,1].set_title('MGE Coverage vs Host Coverage')
    axs[0, 1].set_xlabel('MGE Coverage (log scale)')
    axs[0, 1].set_ylabel('Host Coverage (log scale)')


    # Calculate median cos_sim for each environment and sort
    env_order = df.groupby('environment')['cos_sim'].median().sort_values(ascending=False).index
    sns.boxplot(data=df, x='environment', y='cos_sim', ax=axs[1, 0], hue='environment', legend=False, order=env_order)
    # axs[1,0] set as Cosine Similarity Distribution by Environment
    axs[1, 0].set_xlabel('Environment')
    axs[1, 0].set_ylabel('Cosine Similarity')
    axs[1, 0].tick_params(axis='x', rotation=45, labelsize=9)
    # leave bottom-right empty
    try:
        axs[1, 1].axis('off')
    except Exception:
        pass

    # compute Pearson correlation for coverage on log scale and annotate
    try:
        import numpy as _np
        from scipy.stats import pearsonr
        mask = (df['MGE_cov'] > 0) & (df['host_cov'] > 0)
        if mask.sum() > 1:
            xlog = _np.log10(df.loc[mask, 'MGE_cov'].astype(float))
            ylog = _np.log10(df.loc[mask, 'host_cov'].astype(float))
            r_cov, p_cov = pearsonr(xlog, ylog)
        else:
            r_cov, p_cov = _np.nan, _np.nan
    except Exception:
        try:
            import numpy as _np
            mask = (df['MGE_cov'] > 0) & (df['host_cov'] > 0)
            if mask.sum() > 1:
                xlog = _np.log10(df.loc[mask, 'MGE_cov'].astype(float))
                ylog = _np.log10(df.loc[mask, 'host_cov'].astype(float))
                r_cov = _np.corrcoef(xlog, ylog)[0, 1]
            else:
                r_cov = _np.nan
            p_cov = _np.nan
        except Exception:
            r_cov, p_cov = None, None
    try:
        if r_cov is not None:
            axs[0, 1].text(0.05, 0.95, f"r={r_cov:.2f}\n p={p_cov:.2g}", transform=axs[0, 1].transAxes, va='top', fontsize=9)
    except Exception:
        pass
    # regression lines removed per user request (keep only r/p annotation)

    # improve subplot spacing and save
    try:
        fig.subplots_adjust(wspace=0.35)
    except Exception:
        pass
    plt.tight_layout(pad=2.0)
    plt.savefig(f"{paper_fig_dir}/mge_host_gc_content.pdf")
    plt.close()

def read_drep_cluster(drep_clu_file):
    
    drep_clu_dict = defaultdict(list)
    host_clu_dict = {}
    clu_host_dict = defaultdict(list)
    df = pd.read_csv(drep_clu_file)
    for index, row in df.iterrows():
        # drep_clu_dict[row['genome'][:-3]] = row['secondary_cluster']
        contig = row['genome'][:-3]
        drep_clu_dict[row['secondary_cluster']].append(row['genome'][:-3])
        host_clu_dict[contig] = row['secondary_cluster']
        clu_host_dict[row['secondary_cluster']].append(contig)
    return drep_clu_dict, host_clu_dict, clu_host_dict

def read_mge_cluster(mge_clu_file):
    mge_clu_dict = {}
    clu_mge_dict = defaultdict(list)
    df = pd.read_csv(mge_clu_file, sep="\t", header=None)
    for index, row in df.iterrows():
        cluster = row[0]
        mges = row[1].split(",")
        for mge in mges:
            mge_clu_dict[mge] = cluster
            clu_mge_dict[cluster].append(mge)
    return mge_clu_dict, clu_mge_dict

def count_cross_phylum(whole_G, ctg_taxa_dict, clu_mge_dict=None, rank='phylum'):
    """
    Analyze MGE nodes (virus/plasmid/novel) and count how many are linked to >1
    taxonomic groups at the specified rank (default: 'class').
    """
    cross_count = {'virus': 0, 'plasmid': 0, 'novel': 0}
    type_count = {'virus': 0, 'plasmid': 0, 'novel': 0}
    print ("")
    for node, degree in whole_G.degree:
        node_type = whole_G.nodes[node].get('type')
        if node_type in ['virus', 'plasmid', 'novel']:
            type_count[node_type] += 1
            neighbors = list(whole_G.neighbors(node))
            neighbors_names = []
            groups = set()
            for neighbor in neighbors:
                neigh_type = whole_G.nodes[neighbor].get('type')
                if neigh_type in ['virus', 'plasmid', 'novel']:
                    continue
                neighbor_label = whole_G.nodes[neighbor].get('label', neighbor)
                
                lineage = ctg_taxa_dict.get(neighbor_label, "NA")
                try:
                    group = classify_taxa(lineage, rank)
                except Exception:
                    group = "NA"
                if group and not re.search("Unclassified", str(group)) and group != "NA" and not re.search("unknown", str(group).lower()):
                    groups.add(group)
                    neighbors_names.append((neighbor, group))
            
            if len(groups) > 1:
                clu_info = clu_mge_dict.get(node, 'NA') if clu_mge_dict else 'NA'
                print(f"{node} ({node_type}) is linked to multiple {rank}s: {groups}, mge_cluster: {clu_info}, host contigs: {neighbors_names}")
                cross_count[node_type] += 1

    print(f"Cross-{rank} linked MGEs:")
    print(cross_count)
    print("Number of nodes in each type:")
    print(type_count)

def profile_network(whole_G, ctg_taxa_dict, clu_mge_dict=None):
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
    
    # see which host has linked the most virus (directed: MGE->host, so use predecessors for hosts)
    host_virus_link_count = defaultdict(int)
    for node, degree in whole_G.degree:
        if whole_G.nodes[node]['type'] not in ['virus', 'plasmid', 'novel']:
            # For directed graph: edges are MGE->host, so linked MGEs are predecessors of host
            linked_mges = list(whole_G.predecessors(node))
            for neighbor in linked_mges:
                if whole_G.nodes[neighbor]['type'] == 'virus':
                    host_virus_link_count[node] += 1
    # get the top 10 hosts that linked the most virus
    top_hosts = sorted(host_virus_link_count.items(), key=lambda x: x[1], reverse=True)[:10]
    print("Top 10 hosts linking the most viruses:")
    for host, count in top_hosts:
        node_label = whole_G.nodes[host].get('label', 'Unknown')
        represent_ctg_lineage = ctg_taxa_dict[node_label] if node_label in ctg_taxa_dict else "NA"
        ctg_species = classify_taxa(represent_ctg_lineage, "species")
        ctg_genus = classify_taxa(represent_ctg_lineage, "genus")
        print(f"{host} {node_label}: linked to {count} viruses (Host annotation: {ctg_species}, {ctg_genus})")
    count_cross_phylum(whole_G, ctg_taxa_dict, clu_mge_dict)
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

def analyze_MGEs(gc_df, mge_clu_dict):
    taxon = "s__Enterococcus faecalis"
    target_gc_df = gc_df[gc_df['host_taxa'].str.contains(taxon)]
    print (target_gc_df)
    ## count how many unique MGE clusters, mge_clu_dict
    target_gc_df["MGE_cluster"] = target_gc_df['MGE'].map(mge_clu_dict)
    unique_mge_clusters = target_gc_df['MGE_cluster'].unique()
    print (f"Number of unique MGE clusters linked to {taxon}: {len(unique_mge_clusters)}")
    print (unique_mge_clusters)


def secondary_chr(all_dir):
    genome_dir = "/home/shuaiw/borg/paper/sec_chr/genomes/"
    ### secondary-chromosome
    gc_df = pd.read_csv(f"{paper_fig_dir}/mge_host_gc_cov.csv")
    ## only keep the rows with MGE_type as novel
    target_gc_df = gc_df[gc_df['MGE_type'] == 'novel']
    ## sort by mge_len descending
    target_gc_df = target_gc_df.sort_values(by='mge_len', ascending=False)
    print (target_gc_df.head(20))
    print (len(target_gc_df))
    target_gc_df.to_csv(f"{paper_fig_dir}/secondary_chromosome_gc_cov.csv", index=False)
    for index, row in target_gc_df.iterrows():
        prefix = row['sample']
        mge_obj = My_contig(prefix, all_dir, row['MGE'])
        host_obj = My_contig(prefix, all_dir, row['host'])
        work_dir = genome_dir + "/" +mge_obj.contig
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)
        
        # Create subdirectories if they don't exist
        os.makedirs(f"{work_dir}/MGE_bin", exist_ok=True)
        os.makedirs(f"{work_dir}/host_bin", exist_ok=True)
        os.makedirs(f"{work_dir}/combined_bin", exist_ok=True)
        
        os.system(f"cp {mge_obj.ctg_ref} {work_dir}/MGE_bin")
        os.system(f"cp {host_obj.ctg_ref} {work_dir}/host_bin")
        combined_fasta = f"{work_dir}/combined_bin/{mge_obj.contig}_and_{host_obj.contig}.fa"
        os.system(f"cat {mge_obj.ctg_ref} {host_obj.ctg_ref} > {combined_fasta}")
        print (f"Combined fasta saved to {combined_fasta}")

        ## run checkm2 for each bin folder
        if not os.path.exists(f"{work_dir}/combined_bin/checkm2_output/quality_report.tsv"):
            os.system(f"checkm2 predict --input {work_dir}/MGE_bin --output-directory {work_dir}/MGE_bin/checkm2_output --force --threads 16 -x .fa")
            os.system(f"checkm2 predict --input {work_dir}/host_bin --output-directory {work_dir}/host_bin/checkm2_output --force --threads 16 -x .fa")
            os.system(f"checkm2 predict --input {work_dir}/combined_bin --output-directory {work_dir}/combined_bin/checkm2_output --force --threads 16 -x .fa")

        if  os.path.exists(f"{work_dir}/MGE_bin/checkm2_output/quality_report.tsv"):
            mge_completness, mge_contamination = get_completness(f"{work_dir}/MGE_bin/checkm2_output/quality_report.tsv")
        else:
            mge_completness, mge_contamination = "NA", "NA"
        if os.path.exists(f"{work_dir}/host_bin/checkm2_output/quality_report.tsv"):
             host_completness, host_contamination = get_completness(f"{work_dir}/host_bin/checkm2_output/quality_report.tsv")
        else:
            host_completness, host_contamination = "NA", "NA"
        if os.path.exists(f"{work_dir}/combined_bin/checkm2_output/quality_report.tsv"):
            combined_completness, combined_contamination = get_completness(f"{work_dir}/combined_bin/checkm2_output/quality_report.tsv")
        else:
            combined_completness, combined_contamination = "NA", "NA"


        ## add compleness to df
        target_gc_df.at[index, 'MGE_completness'] = mge_completness
        target_gc_df.at[index, 'Host_completness'] = host_completness
        target_gc_df.at[index, 'Combined_completness'] = combined_completness
        target_gc_df.at[index, 'MGE_contamination'] = mge_contamination
        target_gc_df.at[index, 'Host_contamination'] = host_contamination
        target_gc_df.at[index, 'Combined_contamination'] = combined_contamination
        # break
    print (target_gc_df.head(10))
    target_gc_df.to_csv(f"{paper_fig_dir}/secondary_chromosome_completness.csv", index=False)

def get_completness(checkm_report):
    df = pd.read_csv(checkm_report, sep="\t", header=0)
    for index, row in df.iterrows():
        return row['Completeness'], row["Contamination"]

def component_analysis(connected_components, all_host_clu_lineage_dict, ctg_taxa_dict, paper_fig_dir, clu_host_dict):
    largest_component = connected_components[0]
    # count the nodes number for each type in largest_component
    subgraph = whole_G.subgraph(largest_component)
    type_counts = defaultdict(int)
    for n, d in subgraph.nodes(data=True):
        node_type = d.get('type', 'unknown')
        type_counts[node_type] += 1
    print("Node counts by type in largest connected component:")
    for t, c in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"{t}: {c}")
    
    ## print the degree of the top nodes in largest component
    top_nodes = sorted(subgraph.degree, key=lambda x: x[1], reverse=True)[:10]
    print("Top 10 nodes by degree in largest connected component:")
    for node, degree in top_nodes:
        # also print the type and label of the node
        node_type = subgraph.nodes[node].get('type', 'unknown')
        node_label = subgraph.nodes[node].get('label', node)
        print(f"{node}: {degree} ({node_type}, {node_label})")

    ## according to all_host_clu_lineage_dict, print the species for each host node in largest component
    print("Host nodes and their species in largest connected component:")
    for n, d in subgraph.nodes(data=True):
        node_type = d.get('type', 'unknown')
        if node_type not in ['virus', 'plasmid', 'novel']:
            node_label = d.get('label', n)
            lineage = all_host_clu_lineage_dict.get(n, "NA")
            species = classify_taxa(lineage, "species")
            family = classify_taxa(lineage, "family")
            if isinstance(species, str):
                species = species.replace('s__', '')
            print(f"{n} ({node_label}): {species} (family: {family}) {clu_host_dict[n]}")
    # make a binary matrix for host-MGE linkage in largest component
    host_nodes = [n for n, d in subgraph.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel']]
    mge_nodes = [n for n, d in subgraph.nodes(data=True) if d.get('type') in ['virus', 'plasmid', 'novel']]

    # map host node -> species (use all_host_clu_lineage_dict first, fallback to ctg_taxa_dict via label)
    host_species = {}
    species_counts = defaultdict(int)
    for n in host_nodes:
        label = subgraph.nodes[n].get('label', n)
        lineage = all_host_clu_lineage_dict.get(n, ctg_taxa_dict.get(label, 'NA'))
        try:
            species = classify_taxa(lineage, 'species')
            if not species or re.search('Unclassified', str(species)):
                species = get_detail_taxa_name(lineage)
        except Exception:
            species = get_detail_taxa_name(lineage) if lineage != 'NA' else 'NA'
        host_species[n] = species
        species_counts[species] += 1

    # choose top 14 hosts: prefer hosts from top species and higher degree
    hosts_sorted = sorted(host_nodes, key=lambda n: (species_counts.get(host_species.get(n, 'NA'), 0), subgraph.degree[n]), reverse=True)
    # top_hosts = hosts_sorted[:14]
    top_hosts = hosts_sorted
    skip_heatmap = False
    if not top_hosts:
        print('No host nodes available for heatmap.')
        skip_heatmap = True

    # build human-readable host column labels: original_label (species)
    host_col_labels = []
    for n in top_hosts:
        orig = subgraph.nodes[n].get('label', n)
        sp = host_species.get(n, 'NA')
        host_col_labels.append(f"{orig} ({sp})")

    # map mge node -> original label
    mge_labels = {m: subgraph.nodes[m].get('label', m) for m in mge_nodes}

    # build binary matrix with rows = MGE cluster (node id), cols = Host cluster (node id)
    binary_matrix = pd.DataFrame(0, index=mge_nodes, columns=top_hosts)
    # prepare swapped tick labels: species (original) -> now show species first then original
    host_ticklabels = []
    for n in top_hosts:
        orig = subgraph.nodes[n].get('label', n)
        sp = host_species.get(n, 'NA')
        if isinstance(sp, str):
            sp = sp.replace('s__', '')
        host_ticklabels.append(f"{sp} ({orig})")

    for m in mge_nodes:
        for neigh in subgraph.neighbors(m):
            if neigh in top_hosts:
                binary_matrix.at[m, neigh] = 1

    # move virus rows to the bottom while preserving order otherwise
    non_virus_rows = [r for r in binary_matrix.index if subgraph.nodes[r].get('type') != 'virus']
    virus_rows = [r for r in binary_matrix.index if subgraph.nodes[r].get('type') == 'virus']
    new_index = non_virus_rows + virus_rows
    binary_matrix = binary_matrix.reindex(new_index)

    print('Binary matrix (rows=MGE cluster id, cols=Host cluster id) — showing head:')
    # print(binary_matrix.head(50))
    # build and save heatmap from `binary_matrix`
    if skip_heatmap:
        print('Skipping heatmap generation due to no hosts.')
    else:
        ## store the matrix for potential reuse
        binary_matrix.to_csv(os.path.join(paper_fig_dir, 'largest_component_binary_matrix.csv'))
        figsize = (max(6, len(binary_matrix.columns) * 0.4), max(6, len(binary_matrix.index) * 0.25))
        # attempt hierarchical clustering on non-virus rows, keep viruses at bottom
        try:
            import scipy.cluster.hierarchy as sch
            from scipy.spatial.distance import pdist

            non_virus_rows = [r for r in binary_matrix.index if subgraph.nodes[r].get('type') != 'virus']
            virus_rows = [r for r in binary_matrix.index if subgraph.nodes[r].get('type') == 'virus']

            if len(non_virus_rows) > 1:
                sub_nv = binary_matrix.loc[non_virus_rows]
                row_dist = pdist(sub_nv.values, metric='jaccard') if sub_nv.shape[1] > 0 else pdist(sub_nv.values)
                row_link = sch.linkage(row_dist, method='average')
                row_order = sch.leaves_list(row_link)
            else:
                row_order = list(range(len(non_virus_rows)))

            if len(non_virus_rows) > 0 and binary_matrix.shape[1] > 1:
                col_dist = pdist(binary_matrix.loc[non_virus_rows].values.T, metric='jaccard')
                col_link = sch.linkage(col_dist, method='average')
                col_order = sch.leaves_list(col_link)
            else:
                col_order = list(range(binary_matrix.shape[1]))

            ordered_nonvirus = binary_matrix.loc[non_virus_rows].iloc[row_order, col_order]
            if virus_rows:
                ordered_virus = binary_matrix.loc[virus_rows].iloc[:, col_order]
                binary_ordered = pd.concat([ordered_nonvirus, ordered_virus])
            else:
                binary_ordered = ordered_nonvirus

            ordered_col_ticklabels = [host_ticklabels[i] for i in col_order]

            fig, ax = plt.subplots(figsize=figsize)
            sns.heatmap(binary_ordered, cmap='Greys', cbar=False, linewidths=0.5, linecolor='lightgray',
                        xticklabels=ordered_col_ticklabels, yticklabels=binary_ordered.index.tolist(), ax=ax)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.set_ylabel('MGE cluster')
            ax.set_xlabel('Host cluster')
            ax.set_title('')

            # color MGE row labels by type
            yticklabels = ax.get_yticklabels()
            ylbl_colors = []
            for node in binary_ordered.index:
                t = subgraph.nodes[node].get('type')
                if t == 'plasmid':
                    ylbl_colors.append('blue')
                elif t == 'virus':
                    ylbl_colors.append('red')
                elif t == 'novel':
                    ylbl_colors.append('orange')
                else:
                    ylbl_colors.append('black')

            for lbl_obj, color in zip(yticklabels, ylbl_colors):
                lbl_obj.set_color(color)

            out_heat = os.path.join(paper_fig_dir, 'largest_component_binary_heatmap_from_matrix.pdf')
            plt.tight_layout()
            plt.savefig(out_heat, dpi=300)
            plt.close()
            print(f"Saved clustered heatmap to {out_heat}")
        except Exception:
            # fallback to simple heatmap
            fig, ax = plt.subplots(figsize=figsize)
            sns.heatmap(binary_matrix, cmap='Greys', cbar=False, linewidths=0.5, linecolor='lightgray',
                        xticklabels=host_ticklabels, yticklabels=binary_matrix.index.tolist(), ax=ax)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.set_ylabel('MGE cluster')
            ax.set_xlabel('Host cluster')
            ax.set_title('Binary heatmap: plasmid/virus presence across host hosts')

            yticklabels = ax.get_yticklabels()
            ylbl_colors = []
            for node in binary_matrix.index:
                t = subgraph.nodes[node].get('type')
                if t == 'plasmid':
                    ylbl_colors.append('blue')
                elif t == 'virus':
                    ylbl_colors.append('red')
                elif t == 'novel':
                    ylbl_colors.append('orange')
                else:
                    ylbl_colors.append('black')

            for lbl_obj, color in zip(yticklabels, ylbl_colors):
                lbl_obj.set_color(color)

            out_heat = os.path.join(paper_fig_dir, 'largest_component_binary_heatmap_from_matrix.pdf')
            plt.tight_layout()
            plt.savefig(out_heat, dpi=300)
            plt.close()
            print(f"Saved heatmap to {out_heat} (no clustering)")
        ## Melt the binary matrix, add species and MGE type, and save entries where present==1
        try:
            melt_df = binary_matrix.reset_index().melt(id_vars='index', var_name='host_node', value_name='present')
            melt_df = melt_df.rename(columns={'index': 'MGE_node'})
            # Map labels and types (use mappings available in this scope)
            melt_df['MGE_label'] = melt_df['MGE_node'].map(lambda x: mge_labels.get(x, x))
            melt_df['MGE_type'] = melt_df['MGE_node'].map(lambda x: subgraph.nodes[x].get('type', 'unknown'))
            melt_df['host_label'] = melt_df['host_node'].map(lambda x: subgraph.nodes[x].get('label', x))
            melt_df['host_species'] = melt_df['host_node'].map(lambda x: host_species.get(x, 'NA'))
            # Clean species strings (remove leading 's__' if present)
            melt_df['host_species'] = melt_df['host_species'].astype(str).str.replace('s__', '', regex=False)
            # Keep only positive links (present == 1)
            melt_present = melt_df[melt_df['present'] == 1].copy()
            out_melt = os.path.join(paper_fig_dir, 'largest_component_binary_matrix_melted.csv')
            melt_present.to_csv(out_melt, index=False)
            print(f"Saved melted matrix (present==1) to {out_melt}")
        except Exception as e:
            print(f"Warning: failed to produce melted CSV: {e}")

        ## also save the binary matrix to a csv file
        # binary_matrix.to_csv(os.path.join(paper_fig_dir, 'largest_component_binary_matrix.csv'))

def save_host_list(host_ctg_set, host_list_file):
    with open(host_list_file, "w") as f:
        for ctg in host_ctg_set:
            f.write(f"{ctg}\n")

if __name__ == "__main__":  
    ANI = 99
    meta_file = "/home/shuaiw/MODIFI/assembly_pipe/prefix_table.tab"
    host_list_file =  "../../tmp/figures/multi_env_linkage/linked_host_list.txt"
    mge_clu_file = "/home/shuaiw/borg/paper/MGE/cluster/megablast.cluster.95ani.tsv"
    paper_fig_dir = f"../../tmp/figures/multi_env_linkage/network_{ANI}/"
    # drep_clu_file = f"/home/shuaiw/borg/paper/specificity/dRep_{ANI}_out/data_tables/Cdb.csv"
    drep_clu_file = "/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/dRep_99_out/data_tables/Cdb.csv"
    all_dir = "/home/shuaiw/borg/paper/run2/"
    host_ctg_set = set()
    ## mkdir paper_fig_dir if not exists
    if not os.path.exists(paper_fig_dir):
        os.makedirs(paper_fig_dir)

    
    sample_env_dict = read_metadata(meta_file)
    drep_clu_dict, host_clu_dict, clu_host_dict = read_drep_cluster(drep_clu_file)
    mge_clu_dict, clu_mge_dict = read_mge_cluster(mge_clu_file)
    all_host_clu_lineage_dict = {}
    
    ctg_taxa_dict = get_ctg_taxa(all_dir)
    # """

    whole_G = nx.DiGraph()
    gc_data = []
    cluster_anno_dict = {}
    
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if prefix in [ "ERR5621427_sludge", "ERR5621429_sludge", "ERR5621430_sludge"]:
            continue

        sample_obj = My_sample(prefix, all_dir)
        print(f"Processing {prefix}...")

        MGE_type_dict = sample_obj.read_MGE()
        G, gc_data, cluster_anno_dict, host_clu_lineage_dict = get_edge(cluster_anno_dict, MGE_type_dict, gc_data, 
                                                 sample_env_dict[prefix], prefix, host_clu_dict, 
                                                 mge_clu_dict, sample_obj, ctg_taxa_dict)
        host_ctg_set = collect_host_ctgs(prefix, sample_obj, ctg_taxa_dict,host_ctg_set)
        whole_G = nx.compose(whole_G, G)
        all_host_clu_lineage_dict.update(host_clu_lineage_dict)
    print (f"Total unique host contigs linked to MGEs: {len(host_ctg_set)}")
    save_host_list(host_ctg_set, host_list_file)

    
    plot_network2(whole_G, paper_fig_dir)
    ## assign node colors for GML so Gephi can display them
    assign_node_colors_for_gml(whole_G)
    ## save whole_G to file by gml
    nx.write_gml(whole_G, f"{paper_fig_dir}/whole_network2.gml")
    ## filtered GML: only host nodes in abundant phyla (>5 nodes)
    write_filtered_gml_abundant_phylum(whole_G, f"{paper_fig_dir}/whole_network2_abundant_phylum.gml", min_host_per_phylum=5)
    ## save network to Excel (Nodes + Edges sheets)
    export_network_to_excel(whole_G, f"{paper_fig_dir}/network.xlsx")
    ## save network to JSON (nodes + edges with attributes)
    export_network_to_json(whole_G, f"{paper_fig_dir}/network.json")

    profile_network(whole_G, ctg_taxa_dict, clu_mge_dict)

    gc_df = pd.DataFrame(gc_data, columns=["MGE", "host", "MGE_type", "MGE_gc", "host_gc", 
                                           "cos_sim", "MGE_cov", "host_cov", 
                                           "environment", "sample", "mge_len", "host_taxa", "host_phylum"])
    ## sort gc_df by mge_len descending
    gc_df = gc_df.sort_values(by='MGE_cov', ascending=False)
    ## save gc_df to a csv file
    gc_df.to_csv(f"{paper_fig_dir}/mge_host_gc_cov.csv", index=False)
    plot_gc(gc_df, paper_fig_dir)
    analyze_MGEs(gc_df, mge_clu_dict)

    # secondary_chr(all_dir)
    """

    whole_G = nx.read_gml(f"{paper_fig_dir}/whole_network2.gml")
    print (f"Loaded graph with {whole_G.number_of_nodes()} nodes and {whole_G.number_of_edges()} edges")
    ## priny node number for each type
    type_counts = defaultdict(int)
    for n, d in whole_G.nodes(data=True):
        node_type = d.get('type', 'unknown')
        type_counts[node_type] += 1
    print("Node counts by type in whole graph:")
    for t, c in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"{t}: {c}")
    ## count host nodes number
    host_count = sum(1 for n, d in whole_G.nodes(data=True) if d.get('type') not in ['virus', 'plasmid', 'novel'])
    print(f"Number of host nodes: {host_count}")

    # Output degree of virus and plasmids to CSV (only virus and plasmid nodes)
    export_mge_degrees(whole_G, f"{paper_fig_dir}/virus_plasmid_degrees.csv")
    count_cross_phylum(whole_G, ctg_taxa_dict, clu_mge_dict, rank='phylum')

    ## output the top 5 weakly connected components (directed graph)
    connected_components = sorted(nx.weakly_connected_components(whole_G), key=len, reverse=True)
    for i, component in enumerate(connected_components[:5]):
        print(f"Connected component {i+1} has {len(component)} nodes")
    ## plot the largest connected component


    component_analysis(connected_components, all_host_clu_lineage_dict, ctg_taxa_dict, paper_fig_dir, clu_host_dict)
    """
    # # plot_network2(subgraph, paper_fig_dir, fig_name = "largest_connected_component")
