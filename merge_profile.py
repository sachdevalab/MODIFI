## there are multiple files in the profile dir, we use pandas to read them as a single dataframe
## then we extract the column motif_modified_ratio and merge them into a new dataframe

import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import sys
import argparse
from adjustText import adjust_text
import numpy as np

# sys.setrecursionlimit(20000)

try:
    import fastcluster
except ImportError:
    print("fastcluster not installed. Using scipy for clustering.")

def load_contigs():
    host_file = '/home/shuaiw/Methy/borg_test/borg.csv'
    ## read the contig name in a dict
    contig_dict = {}
    with open(host_file) as f:
        for line in f:
            contig = line.strip()
            contig_dict[contig] = 'borg'
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

def merge_profile(profile_list):
    ## initialize the profiles as df
    profiles = []
    samples = []
    for file in os.listdir(profile_list):
    # for file in profile_list:
        if file.endswith(".csv"):
            profile = pd.read_csv(os.path.join(profile_list, file), sep = ",")
            # profile = pd.read_csv(os.path.join(file), sep = ",")
            ## extract the sample name
            match = re.search(r'(.*?).motifs.profile.csv', file)
            if not match:
                print ("cannot extract sample name from", file)
                continue
            sample_name = match.group(1).split("/")[-1]
            # if sample_name not in borg_contigs:
            #     continue
            # else:
            #     sample_name = sample_name + "_" + borg_contigs[sample_name]
            # print (sample_name)
            samples.append(sample_name)
            # print (sample_name, profile.head())
            motifs_names = profile['motifString']
            # print (motifs_names)
            ## we only want the motif_modified_ratio column
            profile = profile[['motif_modified_ratio']]
            profiles.append(profile['motif_modified_ratio'])
            # if len(samples) > 10:
            #     break
    ## name the columns
    profiles = pd.concat(profiles, axis=1)
    profiles.columns = samples
    ## define the row names as motifs_names
    # print (motifs_names)
    profiles.index = motifs_names
    print (profiles.head())
    ## print the shape of the profiles
    print (profiles.shape)
    
    return profiles

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
        try:
            sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(20, 15))
            plt.savefig(heat_map)
            plt.clf()
        except Exception as e:
            print(f"Failed to create heatmap: {e}")
            plt.figure()
            plt.text(0.5, 0.5, 'Failed to create heatmap', horizontalalignment='center', verticalalignment='center')
            plt.savefig(heat_map)
            plt.clf()

def TSE(df, cluster_fig):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())


    ## reduce dimention using t-SNE
    X_embedded = TSNE(n_components=2).fit_transform(matrix)
    
    clustering = DBSCAN(eps=0.5, min_samples=2).fit(X_embedded)
    # print (clustering.labels_)
    # calculate how many clusters
    n_clusters = len(set(clustering.labels_))
    print (n_clusters, "clusters detected.")
    ## define fig size
    plt.figure(figsize=(10, 10))
    ## plot the cluster result using seaborn
    scatter_plot = sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=clustering.labels_, palette="viridis")

    # ## add labels to each point
    # for i, label in enumerate(df.columns):
    #     scatter_plot.text(X_embedded[i, 0], X_embedded[i, 1], label, fontsize=9)
    ## add labels to each point
    texts = []
    for i, label in enumerate(df.columns):
        texts.append(scatter_plot.text(X_embedded[i, 0], X_embedded[i, 1], label, fontsize=7))
    

    
    ## adjust text to avoid overlap
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))
    plt.savefig(cluster_fig)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get accurate hgt breakpoints", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--all_profiles", type=str, help="<str> separate by space", metavar="\b")
    # required.add_argument("--all_profiles", type=str, help="<str> separate by space", nargs="+", metavar="\b")
    required.add_argument("--heatmap", type=str, help="<str> heatmap file.", metavar="\b")
    required.add_argument("--summary", type=str, help="<str> output motif summary.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())

    # profile_dir = "/home/shuaiw/borg/bench/break2/profiles"
    # profile_dir = sys.argv[1]
    # heat_map = sys.argv[2]
    # total_profile = sys.argv[3]

    profile_list = args['all_profiles']
    heat_map = args['heatmap']
    total_profile = args['summary']
    cluster_fig = "/".join(heat_map.split("/")[:-1]) + "/motif_cluster.pdf"
    # print (profile_dir)
    # profile_dir = "/home/shuaiw/borg/all_test/profiles"
    # heat_map = f"{profile_dir}/../motif_heatmap.pdf"
    # total_profile = f"{profile_dir}/../motif_profile.csv"

    # borg_contigs = load_contigs()
    profiles = merge_profile(profile_list)
    ## save the profiles
    profiles.to_csv(total_profile, index=True)
    # load the profile from the saved file
    # profiles = pd.read_csv("tmp/profiles.csv", index_col=0)
    heatmap(profiles, heat_map)
    TSE(profiles, cluster_fig)


# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/all_test_ccs/profiles         --heatmap  /home/shuaiw/methylation/data/borg/all_test_ccs/motif_heatmap.pdf         --summary /home/shuaiw/methylation/data/borg/all_test_ccs/motif_profile.csv

# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/bench/zymo2/profiles --heatmap  /home/shuaiw/methylation/data/borg/bench/zymo2/motif_heatmap2.pdf  --summary /home/shuaiw/methylation/data/borg/bench/zymo2/motif_profile2.csv
#                                     