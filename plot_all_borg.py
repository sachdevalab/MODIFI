import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import sys

def load_contigs():
    host_file = '/home/shuaiw/Methy/borg_test/borg.csv'
    ## read the contig name in a dict
    contig_dict = {}
    with open(host_file) as f:
        for line in f:
            contig = line.strip()
            contig_dict[contig] = 'borg'
    bin_anno =  "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.bin.csv"
    df = pd.read_csv(bin_anno, header = None)
    for index, row in df.iterrows():
        # print (row)
        if row[1] == "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_UNK":
            continue
        if re.search("Borg", row[1]):
            if row[0] not in contig_dict:
                contig_dict[row[0]] = 'borg'
                print (row[0], 'borg', row[1])
        
    borg_file = '/home/shuaiw/Methy/borg_test/host.csv'
    with open(borg_file) as f:
        for line in f:
            contig = line.strip()
            contig_dict[contig] = 'host'
    control_dict = {}
    control_file = '/home/shuaiw/Methy/borg_test/control.csv'
    with open(control_file) as f:
        for line in f:
            contig = line.strip()
            if contig not in contig_dict:
                control_dict[contig] = 'control'

    contig_dict.update(control_dict)

    return contig_dict


def heatmap(df, heat_map):
    # df = df.T
    # Plot the heatmap with hierarchical clustering
    # sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(30, 60))
    ## check if df is not empty
    if df.empty:
        print ("empty dataframe")
        ## construct an  empty figure
        plt.figure()
        plt.savefig(heat_map)
    else:
        g= sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(10, 6), cbar_kws={'label': 'Intensity'}, yticklabels=1)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=4)  # Y-axis labels
        ## x-axis labels
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=5, rotation=90)
        plt.savefig(heat_map)
        plt.clf()
# datafile = "/home/shuaiw/borg/bench/all_break/motif_profile.csv"
# datafile = "/home/shuaiw/borg/all_test_ccs3/motif_profile.csv"
datafile = "/home/shuaiw/borg/bench/s3/motif_profile.csv"
# datafile = "/home/shuaiw/borg/bench/luis/m84039_230626_221130_s1.hifi_reads.bc2026/motif_profile.csv"
# datafile = "/home/shuaiw/borg/bench/all_subreads/motif_profile2.csv"


# print (name_dict)
df = pd.read_csv(datafile, index_col=0)
## replace the column names with the sample names, and remove the column if not in name_dict
# df = df[[col for col in df.columns if col in name_dict]]
# df.columns = [name_dict[col] for col in df.columns]
# df = df.sample(n=100, axis=0)  #motif number
## remove the rows with all values smaller than 0.1
# print (df.shape)
df = df.T
print (df.head())
contig_dict = load_contigs()
# Filter the DataFrame to keep only the rows with row names in the contig_dict, and rename the rows with the values in the contig_dict
df = df.loc[[row for row in df.index if row in contig_dict]]
df.index = [row + "_" + contig_dict[row] for row in df.index]



## remove the rows with all values smaller than 0.1
df = df.loc[(df > 0.1).any(axis=1)]

## downsample the data
# df = df.sample(n=200, axis=0)  ## contig number
print (df.shape)

## remove the columns with all values smaller than 0.1
df = df.loc[:, (df > 0.1).any(axis=0)]
print (df.shape)

heatmap(df, "tmp/heat_map_borg.pdf")