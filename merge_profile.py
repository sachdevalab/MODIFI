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
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from bin import bin_contigs_to_fastas
from cal_invasion_score import batch_MGE_invade

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

def read_profile_worker(profile_list, file, sample_name,  min_motif_sites=1):
    profile = pd.read_csv(os.path.join(profile_list, file), sep = ",")
    motifs_names = profile['motifString']
    
    profile.loc[profile['motif_modified_num'] < min_motif_sites, 'motif_modified_ratio'] = 0
    ## we only want the motif_modified_ratio column
    profile = profile[['motif_modified_ratio']]
    return profile['motif_modified_ratio'], sample_name, motifs_names

    # profile.loc[profile['motif_modified_num'] < min_motif_sites, 'all_modified_ratio'] = 0
    # ## we only want the motif_modified_ratio column
    # profile = profile[['all_modified_ratio']]
    # return profile['all_modified_ratio'], sample_name, motifs_names

def merge_profile_parallele(profile_list, min_motif_sites=1, threads = 1):
    ## initialize the profiles as df
    profiles = []
    samples = []
    motifs_names = ''
    print ("threads used in profile merge:", threads)
    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for file in os.listdir(profile_list):
        # for file in profile_list:
            if file.endswith(".csv"):
                match = re.search(r'(.*?).motifs.profile.csv', file)
                if not match:
                    print ("cannot extract sample name from", file)
                    continue
                sample_name = match.group(1).split("/")[-1]
                # if sample_name != "E_coli_H10407_1":
                #     continue
                # print (file)
                future = executor.submit(
                    read_profile_worker,
                    profile_list = profile_list,
                    file = file,
                    sample_name = sample_name,
                    min_motif_sites = min_motif_sites,
                )
                futures.append(future)

        for future in tqdm(as_completed(futures), total=len(futures)):
            motif_modified_ratio_column, sample_name, motifs_names = future.result()

            profiles.append(motif_modified_ratio_column)

            samples.append(sample_name)
            motifs_names = motifs_names
    if not profiles:
        return pd.DataFrame()
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

def merge_profile(profile_list, min_motif_sites=1):
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
            ## if motif_modified_num < min_motif_sites, we set the motif_modified_ratio to 0
            ### minimium number of motif sites cutoff, it too less motif sites, no meaning to calculate the ratio
            profile.loc[profile['motif_modified_num'] < min_motif_sites, 'motif_modified_ratio'] = 0
            ## we only want the motif_modified_ratio column
            profile = profile[['motif_modified_ratio']]
            profiles.append(profile['motif_modified_ratio'])

    if not profiles:
        return pd.DataFrame()
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
    bams = [f for f in os.listdir(bam_dir) if f.endswith(".bam") and not f.endswith(".filtered.bam")]
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
    if profile_num > 0:
        ratio_profile_ctg = round(profiles.shape[1]/profile_num,2)
    else:
        ratio_profile_ctg = 0

    f = open(summary_file, "w")
    print (f"No. of contigs with motifs: {profiles.shape[1]}", file=f)
    print (f"Minimum motif methylation fraction: {min_frac}", file=f)
    print (f"Ratio of contigs profiles with motifs: {ratio_profile_ctg}", file=f)
    print (f"No. of motifs: {profiles.shape[0]}", file=f)
    print (f"Number of contigs: {contig_num}", file=f)
    print (f"Number of bams: {bam_num}", file=f)
    print (f"Number of ipds: {ipd_num}", file=f)
    print (f"Number of gffs: {gff_num}", file=f)
    print (f"Number of motif files: {motif_num}", file=f)
    print (f"Number of profiles: {profile_num}", file=f)

    f.close()

def JC_hierarchical_clustering(df, cluster_fig, cutoff=0.45):
    matrix = df.to_numpy()
    ## Trabspose the matrix
    matrix = matrix.T

    matrix = (matrix > 0.5).astype(int) 
    my_linkage = linkage(matrix, method='average', metric='jaccard')
    cluster_labels = fcluster(my_linkage, t=cutoff, criterion='distance')
    n_clusters = len(set(cluster_labels))
    print (n_clusters, "clusters detected in JC.")

    data = []
    for i in range(n_clusters):
        # print ("cluster", i)
        for j in range(len(cluster_labels)):
            if cluster_labels[j] == i:
                # print (df.columns[j])
                data.append([df.columns[j], i])
    cluster_result = pd.DataFrame(data, columns = ['contigs', 'cluster'])
    cluster_result.to_csv(cluster_fig.replace(".pdf", ".j.csv"), index=False)    
    
def heatmap(df, heat_map):
    df = df.T
    # Plot the heatmap with hierarchical clustering
    # sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(30, 60))
    ## check if df is not empty
    if df.empty or df.shape[1] > 200:
        print ("empty or too-large dataframe")
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

def hierarchical_clustering(df, tree_fig, cluster_fig, cutoff=1.6):
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
    try:
        plt.figure(figsize=(20, 15))
        dendrogram(my_linkage, labels=df.columns, orientation='left', leaf_rotation=0)
        # plt.axhline(y=cutoff, color='r', linestyle='--')  # Show the cutoff threshold
        plt.axvline(x=cutoff, color='r', linestyle='--')
        plt.title(f"Dendrogram with Distance Threshold = {cutoff}")
        plt.xlabel("Sample Index")
        plt.ylabel("Cluster Distance")
        plt.savefig(tree_fig)
    except Exception as e:
        print(f"Failed to create dendrogram: {e}")
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

def merge_profile_worker(work_dir, heat_map, profile_list, total_profile, min_frac, whole_ref, plasmid_file, bin_file=None, threads=1, bin_flag = False):
    
    cluster_fig = os.path.join(work_dir, "motif_cluster.pdf")
    pca_fig = os.path.join(work_dir, "motif_pca.pdf")
    tree_fig = os.path.join(work_dir, "motif_tree.pdf")
    summary_file = os.path.join(work_dir, "summary.csv")
    bin_dir = os.path.join(work_dir, "bins", "bin")
    profile_dir = os.path.join(work_dir, "profiles")

    profiles = merge_profile_parallele(profile_list, threads=threads)
    profiles.to_csv(total_profile, index=True)

    # profiles = pd.read_csv(total_profile, index_col=0)
    ## index the prfiles with index and columns
    # if "GATC" in profiles.index:
    #     print(profiles.loc["GATC", "E_coli_H10407_1"])
    # else:
    #     print("Motif 'GATC' not found in profiles.")
    
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    print ("filtered shape 1", profiles.shape)
    # Filter columns where any value is greater than 0.5
    profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]

    print ("filtered shape", profiles.shape)

    if len(profiles) > 0:
        heatmap(profiles, heat_map)

        if len(profiles) > 1 and bin_flag:  ## cluster only if there are more than 1 contig
            # UMAP(profiles, cluster_fig)
            TSE(profiles, cluster_fig)
            # PCA_plot(profiles, pca_fig)
            # hierarchical_clustering(profiles, tree_fig, cluster_fig)
            JC_hierarchical_clustering(profiles, cluster_fig)

            bin_contigs_to_fastas(cluster_fig.replace(".pdf", ".j.csv"), whole_ref, bin_dir)
            
        summary( min_frac, summary_file, profiles)

        # host_dir = os.path.join(work_dir, "hosts")
        # os.makedirs(host_dir, exist_ok = True)
        # if plasmid_file != 'NA' and os.path.exists(plasmid_file):
        #     batch_MGE_invade(plasmid_file, profile_dir, host_dir, bin_file=bin_file, min_frac=0.5, threads=threads)
    else:
        print ("no motif identified")
        ## construct an  empty figure
        plt.figure()
        plt.text(0.5, 0.5, 'No motif identified', horizontalalignment='center', verticalalignment='center')
        plt.savefig(heat_map)
        plt.clf()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get accurate hgt breakpoints", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--all_profiles", type=str, help="<str> profile dir", metavar="\b")
    # required.add_argument("--all_profiles", type=str, help="<str> separate by space", nargs="+", metavar="\b")
    required.add_argument("--heatmap", type=str, help="<str> heatmap file.", metavar="\b")
    required.add_argument("--summary", type=str, help="<str> output motif summary.", metavar="\b")
    required.add_argument("--whole_ref", type=str, help="<str> raw ref with contigs.", metavar="\b")
    required.add_argument("--min_ctg_len", type=int, default = 1000, help="minimal contig length.", metavar="\b")
    optional.add_argument("--min_frac", type=float, default=0.5, help="<float> minimum fraction of methylation to keep the motif.")
    optional.add_argument("--plasmid_file", type=str, help="<str> *_plasmid_summary.tsv by genomad.", default = 'NA', metavar="\b")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())

    min_frac = args['min_frac']
    profile_list = args['all_profiles']
    heat_map = args['heatmap']
    total_profile = args['summary']
    plasmid_file = args['plasmid_file']
    whole_ref = args['whole_ref']

    work_dir = "/".join(heat_map.split("/")[:-1])
    merge_profile_worker(work_dir, heat_map, profile_list, total_profile, min_frac, whole_ref, plasmid_file)







# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/all_test_ccs/profiles         --heatmap  /home/shuaiw/methylation/data/borg/all_test_ccs/motif_heatmap.pdf         --summary /home/shuaiw/methylation/data/borg/all_test_ccs/motif_profile.csv

# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/bench/zymo2/profiles --heatmap  /home/shuaiw/methylation/data/borg/bench/zymo2/motif_heatmap2.pdf  --summary /home/shuaiw/methylation/data/borg/bench/zymo2/motif_profile2.csv
# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/bench/all_break/profiles --heatmap /home/shuaiw/methylation/data/borg/bench/all_break/test.pdf  --summary /home/shuaiw/methylation/data/borg/bench/all_break/motif_profile.csv    
# python /home/shuaiw/Methy/merge_profile.py --all_profiles /home/shuaiw/methylation/data/borg/bench/zymo6_NM200/profiles --heatmap  /home/shuaiw/methylation/data/borg/bench/zymo6_NM200/motif_heatmap2.pdf  --summary /home/shuaiw/methylation/data/borg/bench/zymo6_NM200/motif_profile2.csv                                