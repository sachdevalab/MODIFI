## there are multiple files in the profile dir, we use pandas to read them as a single dataframe
## then we extract the column motif_modified_ratio and merge them into a new dataframe

import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster

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
            samples.append(sample_name)
            # print (sample_name, profile.head())
            motifs_names = profile['motifString']
            # print (motifs_names)
            ## we only want the motif_modified_ratio column
            profile = profile[['motif_modified_ratio']]
            profiles.append(profile['motif_modified_ratio'])

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

def summary( min_frac, summary_file, profiles):
    
    work_dir = "/".join(summary_file.split("/")[:-1])
    ## count number of bams
    bam_dir = os.path.join(work_dir, "bams")
    ## count the number of bams
    bams = [f for f in os.listdir(bam_dir) if f.endswith(".bam")]
    bam_num = len(bams)
    ## count the number of contigs
    contig_file = os.path.join(work_dir, "contigs")
    contigs = [f for f in os.listdir(contig_file) if f.endswith(".fa")]
    contig_num = len(contigs)
    ## count the number of ipds
    ipd_file = os.path.join(work_dir, "ipd")
    ipds = [f for f in os.listdir(ipd_file) if f.endswith(".ipd1.csv")]
    ipd_num = len(ipds)
    ## count the number of profiles
    profile_file = os.path.join(work_dir, "profiles")
    profiles_list = [f for f in os.listdir(profile_file) if f.endswith(".motifs.profile.csv")]
    profile_num = len(profiles_list)
    ## count the number of motif files
    motif_file = os.path.join(work_dir, "motifs")
    motifs = [f for f in os.listdir(motif_file) if f.endswith(".motifs.csv")]
    motif_num = len(motifs)
    ## count the number of gff files
    gff_file = os.path.join(work_dir, "gffs")
    gffs = [f for f in os.listdir(gff_file) if f.endswith(".gff") and not f.endswith(".reprocess.gff")]
    gff_num = len(gffs)

    f = open(summary_file, "w")
    print (f"No. of contigs with motifs: {profiles.shape[1]}", file=f)
    print (f"Minimum motif methylation fraction: {min_frac}", file=f)
    print (f"Ratio of contigs profiles with motifs: {round(profiles.shape[1]/profile_num,2)}", file=f)
    print (f"No. of motifs: {profiles.shape[0]}", file=f)
    print (f"Number of contigs: {contig_num}", file=f)
    print (f"Number of bams: {bam_num}", file=f)
    print (f"Number of ipds: {ipd_num}", file=f)
    print (f"Number of gffs: {gff_num}", file=f)
    print (f"Number of motif files: {motif_num}", file=f)
    print (f"Number of profiles: {profile_num}", file=f)

    f.close()

    
    


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

def hierarchical_clustering(df, tree_fig, cutoff=1.6):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())

    my_linkage = linkage(matrix, method='average', metric='euclidean')
    cluster_labels = fcluster(my_linkage, t=cutoff, criterion='distance')

    ## calculate how many clusters
    n_clusters = len(set(cluster_labels))
    print (n_clusters, "clusters detected hierarchical_clustering.")
    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(cluster_labels)):
            if cluster_labels[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(cluster_fig.replace(".pdf", ".h.csv"), index=False)
    # Plot dendrogram
    plt.figure(figsize=(20, 15))
    dendrogram(my_linkage, labels=df.columns, orientation='left', leaf_rotation=0)
    plt.axhline(y=cutoff, color='r', linestyle='--')  # Show the cutoff threshold
    plt.title(f"Dendrogram with Distance Threshold = {cutoff}")
    plt.xlabel("Sample Index")
    plt.ylabel("Cluster Distance")
    plt.savefig(tree_fig)
    print ("hierarchical clustering done.")

def UMAP(df, cluster_fig):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())
    # print (matrix)

    try:
        import umap
        X_embedded = umap.UMAP().fit_transform(matrix, 
                                               n_neighbors=1,
                                               min_dist=0.1,
                                               metric='euclidean')
        # import umap.plot
        # umap.plot.points(X_embedded)
        # ## save the umap result
        # plt.savefig("/tmp/umap.pdf")
    except Exception as e:
        print(f"Failed to create UMAP: {e}")
        return
    
    clustering = DBSCAN(eps=0.4, min_samples=1).fit(X_embedded)
    # print (clustering.labels_)
    # calculate how many clusters
    n_clusters = len(set(clustering.labels_))
    print (n_clusters, "clusters detected in UMAP.")

    ## output the cluster result, the elements with same cluster label are output together
    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(clustering.labels_)):
            if clustering.labels_[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    ## save the cluster result
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(cluster_fig.replace(".pdf", ".u.csv"), index=False)

def TSE(df, cluster_fig):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())
    # print (matrix)

    try:
        ## reduce dimention using t-SNE
        X_embedded = TSNE(n_components=2).fit_transform(matrix)
    except Exception as e:
        print(f"Failed to create t-SNE: {e}")
        return
    
    clustering = DBSCAN(eps=0.2, min_samples=1).fit(X_embedded)
    # print (clustering.labels_)
    # calculate how many clusters
    n_clusters = len(set(clustering.labels_))
    print (n_clusters, "clusters detected in TSNE.")

    ## output the cluster result, the elements with same cluster label are output together
    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(clustering.labels_)):
            if clustering.labels_[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    ## save the cluster result
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(cluster_fig.replace(".pdf", ".csv"), index=False)
    ## define fig size
    plt.figure(figsize=(10, 10))
    ## plot the cluster result using seaborn
    scatter_plot = sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=clustering.labels_, palette="viridis")
    
    ## adjust text to avoid overlap
    ## if df column number is too large, the adjust_text function may not work well, skip it
    if len(df.columns) < 100:
        texts = []
        for i, label in enumerate(df.columns):
            texts.append(scatter_plot.text(X_embedded[i, 0], X_embedded[i, 1], label, fontsize=7))
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))
    plt.savefig(cluster_fig)

def PCA_plot(df, pca_fig):
    print ("start PCA...")
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    ## transpose the profiles
    # profiles = profiles.T
    ## add pseudo values to zero values
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    ## zero values are set to small random pseudovalues in the (−0.2, +0.2)
    mask = matrix == 0
    matrix[mask] = np.random.uniform(-0.2, 0.2, mask.sum())
    X_embedded = pca.fit_transform(matrix)

    ## cluster the pca result
    clustering = DBSCAN(eps=0.1, min_samples=1).fit(X_embedded)
    # print (clustering.labels_)
    ## save the cluster result in a dataframe, and plot it like in tse function
    n_clusters = len(set(clustering.labels_))
    print (n_clusters, "clusters detected after PCA.")
    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(clustering.labels_)):
            if clustering.labels_[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    ## save the cluster result
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(cluster_fig.replace(".pdf", ".p.csv"), index=False)

    scatter_plot = sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=clustering.labels_, palette="viridis")
    
    ## if df column number is too large, the adjust_text function may not work well, skip it
    if len(df.columns) < 100:
        texts = []
        for i, label in enumerate(df.columns):
            # print (i, label)
            if i < len(X_embedded):
                texts.append(scatter_plot.text(X_embedded[i, 0], X_embedded[i, 1], label, fontsize=7))
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray'))
    plt.savefig(pca_fig)


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

    min_frac = 0.5

    profile_list = args['all_profiles']
    heat_map = args['heatmap']
    total_profile = args['summary']
    cluster_fig = "/".join(heat_map.split("/")[:-1]) + "/motif_cluster.pdf"
    pca_fig = "/".join(heat_map.split("/")[:-1]) + "/motif_pca.pdf"
    tree_fig = "/".join(heat_map.split("/")[:-1]) + "/motif_tree.pdf"
    summary_file = "/".join(heat_map.split("/")[:-1]) + "/summary.csv"

    profiles = merge_profile(profile_list)
    ## save the profiles
    profiles.to_csv(total_profile, index=True)
    
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    # Filter columns where any value is greater than 0.5
    profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]

    print (profiles.shape)

    if len(profiles) > 0:
        heatmap(profiles, heat_map)
        UMAP(profiles, cluster_fig)
        TSE(profiles, cluster_fig)
        PCA_plot(profiles, pca_fig)
        hierarchical_clustering(profiles, tree_fig)
    else:
        print ("no motif identified")
        ## construct an  empty figure
        plt.figure()
        plt.text(0.5, 0.5, 'No motif identified', horizontalalignment='center', verticalalignment='center')
        plt.savefig(heat_map)
        plt.clf()

    summary( min_frac, summary_file, profiles)


# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/all_test_ccs/profiles         --heatmap  /home/shuaiw/methylation/data/borg/all_test_ccs/motif_heatmap.pdf         --summary /home/shuaiw/methylation/data/borg/all_test_ccs/motif_profile.csv

# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/bench/zymo2/profiles --heatmap  /home/shuaiw/methylation/data/borg/bench/zymo2/motif_heatmap2.pdf  --summary /home/shuaiw/methylation/data/borg/bench/zymo2/motif_profile2.csv
# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/bench/all_break/profiles --heatmap /home/shuaiw/methylation/data/borg/bench/all_break/test.pdf  --summary /home/shuaiw/methylation/data/borg/bench/all_break/motif_profile.csv                                    