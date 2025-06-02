
import pandas as pd
import networkx as nx
import re
import os
import matplotlib.pyplot as plt



    # Read the host summary file
    


host_sum_file = "/home/shuaiw/methylation/data/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/host_summary.csv"
df = pd.read_csv(host_sum_file)

## keep the rows with final_score > 0.5
df = df[df['final_score'] > 0.5]
## plot dot plot using seaborn, x is MGE_gc, y is host_gc
plt.figure(figsize=(6, 6))
plt.scatter(df['MGE_gc'], df['host_gc'], alpha=0.5)
plt.title('MGE GC content vs Host GC content')
plt.xlabel('MGE GC content')
plt.ylabel('Host GC content')
plt.grid(True)
plt.savefig('../../tmp/results/mge_host_gc_content.png')
plt.close()

## plot the distribution of cos_sim
plt.figure(figsize=(6, 6))
plt.hist(df['cos_sim'], bins=50, alpha=0.7, color='blue')
plt.title('Distribution of Cosine Similarity')
plt.xlabel('Cosine Similarity')
plt.ylabel('Frequency')
plt.grid(True)
plt.savefig('../../tmp/results/cos_sim_distribution.png')
plt.close()