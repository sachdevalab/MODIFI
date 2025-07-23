import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def heatmap(df, heat_map, plasmid_list, host_list):
    df = df.T

    # Assign type to each sequence
    type_labels = []
    for seq in df.index:
        if seq in plasmid_list:
            type_labels.append('plasmid')
        elif seq in host_list:
            type_labels.append('host')
        else:
            type_labels.append('host')

    # Map types to colors
    type_palette = {'plasmid': '#1f77b4', 'host': '#ff7f0e', 'unknown': '#bbbbbb'}
    row_colors = pd.Series(type_labels, index=df.index).map(type_palette)

    # Plot clustered heatmap
    g = sns.clustermap(
        df,
        row_colors=row_colors,
        method='average',
        metric='euclidean',
        cmap='viridis',
        figsize=(25, 15)
    )
    plt.savefig(heat_map)
    plt.close()


def get_new_host(plasmid_list_file):
    plasmid_host_dict = {}
    plasmid_list = []
    host_list = []
    df = pd.read_csv(plasmid_list_file)
    for index, row in df.iterrows():
        plasmid = row['seq_name']
        host_str = row['host']
        plasmid_list.append(plasmid)
        host_list += host_str.split(";")

    host_list = list(set(host_list))
    return plasmid_list, host_list


motif_profile = "/home/shuaiw/borg/paper/linkage/m64004_210929_143746.p100/motif_profile.csv"
plasmid_list_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list"
heat_map = "./motif_heatmap.pdf"


df = pd.read_csv(motif_profile, index_col=0)
plasmid_list, host_list = get_new_host(plasmid_list_file)
heatmap(df, heat_map, plasmid_list, host_list)