import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import sys



def heatmap(df, heat_map):
    df = df.T
    # Plot the heatmap with hierarchical clustering
    # sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(30, 60))
    ## check if df is not empty
    if df.empty:
        print ("empty dataframe")
        ## construct an  empty figure
        plt.figure()
        plt.savefig(heat_map)
    else:
        sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(15, 8), cbar_kws={'label': 'Intensity'})
        plt.savefig(heat_map)
        plt.clf()
datafile = "/home/shuaiw/methylation/data/borg/bench/mock/test.profile"
ref = "/home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa"

## use biopython to read the reference
from Bio import SeqIO
name_dict = {}
for record in SeqIO.parse(ref, "fasta"):
    contig = record.id
    contig_len = len(record.seq)
    if contig_len < 100000:
        continue
    ## get the annotation from the contig
    match = re.search(r'Mixed sample (.*?),', record.description)
    if not match:
        print ("cannot extract sample name from", record.description)
        continue
    name = match.group(1)
    name_dict[contig] = name


# print (name_dict)
df = pd.read_csv(datafile, index_col=0)
## replace the column names with the sample names, and remove the column if not in name_dict
df = df[[col for col in df.columns if col in name_dict]]
df.columns = [name_dict[col] for col in df.columns]
df = df.T
print (df.head())

heatmap(df, "tmp/heat_map.pdf")