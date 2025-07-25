import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
# ...existing code...


def r_output(df):
    data= []
    for index, row in df.iterrows():
        
        seq = row.name
        for col in df.columns:
            value = row[col]
            ## if cannot convert to float, skip
            print (row)
            try:
                value = float(value)
            except ValueError:
                print (row)
                break
                
            # print (value)
            # if col == "GTAC":
            #     print (f"###{seq}&&\t{col}&&\t{value}>>>")
            data.append([seq, col, value])
    r_df = pd.DataFrame(data, columns=['seq', 'motif', 'value'])
    r_df.to_csv("heatmap.csv", index=False)

def heatmap(df, heat_map, plasmid_list, host_list):
    r_output(df)
    df = df.T
    
    ## output the df index to list
    seq_list = sorted(df.index.tolist())
    ### sort df by seq_list
    df = df.reindex(seq_list)
    # print (df)
    ## set the value to 0 if the value is less than 0.5
    # df[df < 0.4] = 0

    # r_output(df)


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

    # Create custom legend patches
    handles = [
        mpatches.Patch(color='#1f77b4', label='plasmid'),
        mpatches.Patch(color='#ff7f0e', label='host'),
        # mpatches.Patch(color='#bbbbbb', label='unknown')
    ]
    # Plot clustered heatmap
    g = sns.clustermap(
        df,
        row_colors=row_colors,
        method='average',
        metric='euclidean',
        cmap='viridis',
        figsize=(25, 15),
        row_cluster=False,  # Only cluster columns
        col_cluster=True
    )

    # Add the legend to the right of the heatmap
    g.ax_heatmap.legend(
        handles=handles,
        title="Type",
        loc='upper right',
        borderaxespad=0.
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