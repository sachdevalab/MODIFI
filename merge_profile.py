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

def merge_profile(profile_dir):
    
    ## initialize the profiles as df
    profiles = []
    samples = []
    for file in os.listdir(profile_dir):
        if file.endswith(".csv"):
            profile = pd.read_csv(os.path.join(profile_dir, file), sep = ",")
            ## extract the sample name
            match = re.search(r'(.*?).motifs.profile.csv', file)
            if not match:
                print ("cannot extract sample name from", file)
                continue
            sample_name = match.group(1)
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
        sns.clustermap(df, method='average', metric='euclidean', cmap='viridis')
        plt.savefig(heat_map)
        plt.clf()

def TSE(df):
    matrix = df.to_numpy()
    print (matrix.shape)
    ## reduce dimention using t-SNE
    
    X_embedded = TSNE(n_components=2).fit_transform(matrix)
    
    clustering = DBSCAN(eps=0.5, min_samples=2).fit(X_embedded)
    # print (clustering.labels_)
    ## define fig size
    plt.figure(figsize=(10, 10))
    ## plot the cluster result using seaborn
    sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=clustering.labels_)
    # plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c = clustering.labels_)
    ## save the figure in pdf
    plt.savefig("tmp/tse.pdf")

if __name__ == "__main__":
    # profile_dir = "/home/shuaiw/borg/bench/break2/profiles"
    profile_dir = sys.argv[1]
    # profile_dir = "/home/shuaiw/borg/all_test/profiles"
    # heat_map = f"{profile_dir}/../motif_heatmap.pdf"
    # total_profile = f"{profile_dir}/../motif_profile.csv"
    heat_map = sys.argv[2]
    total_profile = sys.argv[3]
    # borg_contigs = load_contigs()
    profiles = merge_profile(profile_dir)
    ## save the profiles
    profiles.to_csv(total_profile, index=True)
    # load the profile from the saved file
    # profiles = pd.read_csv("tmp/profiles.csv", index_col=0)
    heatmap(profiles, heat_map)
    # TSE(profiles)