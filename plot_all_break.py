import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import sys



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
        ## define the figure size
        plt.figure(figsize=(30, 60))
        g= sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(10, 12), cbar_kws={'label': 'Intensity'}, yticklabels=1)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=4)  # Y-axis labels
        ## x-axis labels
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=5, rotation=90)
        plt.savefig(heat_map)
        plt.clf()
# datafile = "/home/shuaiw/borg/bench/all_break/motif_profile.csv"
datafile = "/home/shuaiw/borg/all_test_ccs2/motif_profile.csv"


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


first_column = ['SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_14889_L', 'SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_33175_L']
# fai = "/home/shuaiw/methylation/data/borg/contigs/all_break.contigs.fa.fai"
fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
with open(fai) as f:
    for line in f:
        first_column.append(line.split()[0])
# first_column = first_column[:500]

# Filter the DataFrame to keep only the rows with row names in the contig names list
df = df[df.index.isin(first_column)]
## remove the rows with all values smaller than 0.1
df = df.loc[(df > 0.5).any(axis=1)]

## downsample the data
# df = df.sample(n=200, axis=0)  ## contig number
print (df.shape)

## remove the columns with all values smaller than 0.1
df = df.loc[:, (df > 0.5).any(axis=0)]
print (df.shape)

heatmap(df, "tmp/heat_map_all_ccs.pdf")