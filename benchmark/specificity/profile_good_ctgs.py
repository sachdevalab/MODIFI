import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics.pairwise import cosine_similarity
import re

from Bio.Seq import Seq

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_unique_motifs


def get_best_ctg(depth_file, fai, min_len = 1000000):
    """
    Get the best contig based on length from a fasta file.
    """
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    length_dict = {}
    for index, row in depth_df.iterrows():

        if row['depth'] >= 10:
            good_depth[row['contig']] = row['depth']
    best_ctgs = {}
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            length_dict[ctg] = length
            if ctg[-1] == "C" and length >= min_len and ctg in good_depth:
                best_ctgs[ctg] = length
    # print (f"Total {len(best_ctgs)} contigs with length >= {min_len} found.")
    return good_depth, length_dict

def get_drep_motifs():
    drep_motifs = "/home/shuaiw/borg/bench/soil/run1/test.motifs_drep.csv"
    drep_motif_df = pd.read_csv(drep_motifs)
    drep_motifs_dict = {}
    for index, row in drep_motif_df.iterrows():
        motif = row['motifString']
        if motif not in drep_motifs_dict:
            drep_motifs_dict[motif] = 1
    return drep_motifs_dict

def plot_strain_heatmap(profile_file, min_frac=0.4):
    profiles = pd.read_csv(profile_file, index_col=0)
    ## only keep the columns that are in the best contigs
    best_ctgs = get_best_ctg()
    multiple_strain_ctg, ctg_phylum = read_drep_cluster()
    ## retain only the contigs that are in multiple strain clusters
    best_ctgs = [ctg for ctg in best_ctgs if ctg in multiple_strain_ctg]
    print (f"Total {len(best_ctgs)} contigs")
    
    profiles = profiles.loc[:, profiles.columns.isin(best_ctgs)]
    print ("original shape", profiles.shape)                                    
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    print ("filtered shape 1", profiles.shape)
    print (profiles)
    drep_motifs_dict = get_drep_motifs()


    ## filter profiles to keep only the rows with motifString  that are in drep_motifs_dict
    profiles = profiles[profiles.index.isin(drep_motifs_dict.keys())]
    print ("filtered shape 2", profiles.shape)


    df = profiles.T
    # ctg_phylum = get_phylum(best_ctgs)

    phylum_colors = sns.color_palette("husl", len(set(ctg_phylum.values())))
    phylum_color_dict = {phylum: color for phylum, color in zip(set(ctg_phylum.values()), phylum_colors)}
    phylum_colors_list = [phylum_color_dict[ctg_phylum[ctg]] for ctg in df.index]

    # Use clustermap instead of heatmap
    g = sns.clustermap(
        df,
        method='average',
        metric='euclidean',
        cmap='viridis',
        figsize=(30, 15),
        row_colors=phylum_colors_list,
        xticklabels=True,
        yticklabels=True
    )

    # Add legend for phylum colors
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=phylum_color_dict[phylum], label=phylum) for phylum in phylum_color_dict]
    plt.legend(handles=handles, title="Strain", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


    plt.savefig("../../tmp/results/heatmap_strain_cluster.png", dpi=300, bbox_inches='tight')
    plt.clf()

def plot_heatmap(profile_file, min_frac=0.4):
    profiles = pd.read_csv(profile_file, index_col=0)
    ## only keep the columns that are in the best contigs
    best_ctgs = get_best_ctg()
    
    profiles = profiles.loc[:, profiles.columns.isin(best_ctgs)]
    print ("original shape", profiles.shape)                                    
    profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    print ("filtered shape 1", profiles.shape)
    print (profiles)
    drep_motifs_dict = get_drep_motifs()


    ## filter profiles to keep only the rows with motifString  that are in drep_motifs_dict
    profiles = profiles[profiles.index.isin(drep_motifs_dict.keys())]
    print ("filtered shape 2", profiles.shape)

    # Filter columns where any value is greater than 0.5
    # profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]
    # print ("filtered shape", profiles.shape)

    df = profiles.T
    ctg_phylum = get_phylum(best_ctgs)

    phylum_colors = sns.color_palette("husl", len(set(ctg_phylum.values())))
    phylum_color_dict = {phylum: color for phylum, color in zip(set(ctg_phylum.values()), phylum_colors)}
    phylum_colors_list = [phylum_color_dict[ctg_phylum[ctg]] for ctg in df.index]

    # Use clustermap instead of heatmap
    g = sns.clustermap(
        df,
        method='average',
        metric='euclidean',
        cmap='viridis',
        figsize=(30, 15),
        row_colors=phylum_colors_list,
        xticklabels=True,
        yticklabels=True
    )

    # Add legend for phylum colors
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=phylum_color_dict[phylum], label=phylum) for phylum in phylum_color_dict]
    plt.legend(handles=handles, title="Phylum", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


    plt.savefig("../../tmp/results/heatmap.png", dpi=300, bbox_inches='tight')
    plt.clf()

    # plot_drep_similarity_vs_cos_sim(df)

def plot_drep_similarity_vs_cos_sim(df):
    ## compute pairwise cosine similarity matrix
    cosine_sim_matrix = pd.DataFrame(
        cosine_similarity(df),
        index=df.index,
        columns=df.index
    )
    drep_sim_dict = get_ctg_sim()
    ## for each pair of contigs, get the Euclidean distance and dRep similarity
    drep_sim_list = []
    for ctg1 in df.index:
        for ctg2 in df.index:
            if ctg1 != ctg2:
                euclidean_dist = cosine_sim_matrix.loc[ctg1, ctg2]
                if ctg1 + "&" + ctg2 not in drep_sim_dict:
                    print (f"Not Found dRep similarity for {ctg1} and {ctg2}")
                    continue
                drep_sim = drep_sim_dict[ctg1 + "&" + ctg2]
                drep_sim_list.append((ctg1, ctg2, euclidean_dist, drep_sim))
    drep_sim_df = pd.DataFrame(drep_sim_list, columns=['ctg1', 'ctg2', 'cos_sim', 'drep_sim'])
    ### plot scatter plot of dRep similarity vs Euclidean distance
    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=drep_sim_df, x='cos_sim', y='drep_sim', alpha=0.5)
    plt.xlabel('Cosine Similarity')
    plt.ylabel('dRep Similarity')
    plt.savefig("../../tmp/results/drep_similarity_vs_cos_sim.png", dpi=300, bbox_inches='tight')
    # # print (drep_sim_dict)

def get_ctg_sim():
    drep_sim_dict = {}
    drep_sim = "/home/shuaiw/borg/contigs/dRep/dRep_out4/data_tables/Mdb.csv"
    df = pd.read_csv(drep_sim)
    for index, row in df.iterrows():
        ctg1 = row['genome1'][:-6]
        ctg2 = row['genome2'][:-6]
        drep_sim_dict[ctg1 + "&" + ctg2] = row['similarity']
        drep_sim_dict[ctg2 + "&" + ctg1] = row['similarity']
    return drep_sim_dict

def get_ctg_sim2():
    drep_sim_dict = {}
    drep_sim = "/home/shuaiw/borg/contigs/fastANI.output"
    df = pd.read_csv(drep_sim, header =None, sep="\t")
    for index, row in df.iterrows():
        ## get first colunm as ctg1 and second column as ctg2

        ctg1 = row[0][:-6].split("/")[-1]
        ctg2 = row[1][:-6].split("/")[-1]
        drep_sim_dict[ctg1 + "&" + ctg2] = row[2]
        drep_sim_dict[ctg2 + "&" + ctg1] = row[2]
    # print (drep_sim_dict)
    return drep_sim_dict

def read_drep_cluster(drep_clu_file, best_ctgs):
    
    drep_clu_dict = {}
    df = pd.read_csv(drep_clu_file)
    for index, row in df.iterrows():
        drep_clu_dict[row['genome'][:-6]] = row['secondary_cluster']
    ## count the number of contigs in each cluster
    clu_count = {}
    for ctg, clu in drep_clu_dict.items():
        if ctg not in best_ctgs:
            continue
        if clu not in clu_count:
            clu_count[clu] = 0
        clu_count[clu] += 1
    multiple_strain_ctg = {}
    ctg_phylum = {}
    clu_set = set()
    for ctg, clu in drep_clu_dict.items():
        if clu not in clu_count:
            continue
        if clu_count[clu] > 1:
            multiple_strain_ctg[ctg] = clu
            ctg_phylum[ctg] = clu
            clu_set.add(clu)
    print ("multiple_strain_ctg", multiple_strain_ctg)
    return multiple_strain_ctg, ctg_phylum, len(clu_set)

def jaccard():
    # drep_clu_file = "/home/shuaiw/methylation/data/borg/contigs/dRep/dRep_out/data_tables/Cdb.csv"
    # profile_file = "/home/shuaiw/borg/bench/soil/run1/motif_profile.csv"
    # depth_file = "/home/shuaiw/borg/bench/soil/run1/mean_depth.csv"

    prefix = "ERR12723529"
    work_dir = "/home/shuaiw/borg/paper/run"
    drep_clu_file = os.path.join(work_dir, prefix, "dRep_out", "data_tables", "Cdb.csv")
    profile_file = os.path.join(work_dir, prefix, f"{prefix}_methylation", "motif_profile.csv")
    depth_file = os.path.join(work_dir, prefix, f"{prefix}_methylation", "mean_depth.csv")

    multiple_strain_ctg, ctg_phylum = read_drep_cluster(drep_clu_file)
    ## read motif profile file
    profiles = pd.read_csv(profile_file, index_col=0)
    ## binary it by 0.5
    profiles_bin = (profiles >= 0.5).astype(int)
    profiles_bin = profiles_bin.T
    ## calculate Jaccard similarity in the same cluster and different clusters
    jaccard_similarities = []
    ctg_list = list(multiple_strain_ctg.keys())
    # print (profiles_bin)
    for i in range(len(ctg_list)):
        for j in range(i + 1, len(ctg_list)):
            ctg1 = ctg_list[i]
            ctg2 = ctg_list[j]
            ## skip if ctg1 not in profiles_bin or ctg2 not in profiles_bin
            if ctg1 not in profiles_bin.index or ctg2 not in profiles_bin.index:
                continue
            if ctg1 == ctg2:
                continue
            # print (profiles_bin.loc[ctg1])
            intersection = np.sum((profiles_bin.loc[ctg1] & profiles_bin.loc[ctg2]).astype(int))
            union = np.sum((profiles_bin.loc[ctg1] | profiles_bin.loc[ctg2]).astype(int))
            if union == 0:
                continue
            jaccard_similarity = intersection / union if union > 0 else 0
            if multiple_strain_ctg[ctg1] != multiple_strain_ctg[ctg2]:
                # Calculate Jaccard similarity
                jaccard_similarities.append([jaccard_similarity, "diff"])
            else:
                jaccard_similarities.append([jaccard_similarity, "same"])
    print (f"Total {len(jaccard_similarities)} Jaccard similarities calculated.")
    ## plot the Jaccard similarities
    df_jaccard = pd.DataFrame(jaccard_similarities, columns=['jaccard_similarity', 'cluster_type'])
    plt.figure(figsize=(7, 6))
    sns.boxplot(data=df_jaccard, x='cluster_type', y='jaccard_similarity', palette='Set2')
    plt.title('Jaccard Similarity of Motifs in Multiple Strain Clusters')
    plt.xlabel('Cluster Type')
    plt.ylabel('Jaccard Similarity')
    plt.savefig("../../tmp/results/jaccard_similarity.png", dpi=300, bbox_inches='tight')
    plt.clf()   

def jaccard_batch():
    all_dir = "/home/shuaiw/borg/paper/run2/"
    jaccard_similarities = []
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        print (f"Processing {prefix}...")
        # work_dir = f"{all_dir}/{prefix}/{prefix}_methylation3"

        # prefix = "ERR12723529"
        work_dir = "/home/shuaiw/borg/paper/run2"
        # drep_clu_file = os.path.join(work_dir, prefix, "dRep_out_99", "data_tables", "Cdb.csv")
        # drep_clu_file = os.path.join(work_dir, prefix, "dRep_out", "data_tables", "Cdb.csv")
        drep_clu_file = os.path.join(work_dir, prefix, "dRep_out_97", "data_tables", "Cdb.csv")
        profile_file = os.path.join(work_dir, prefix, f"{prefix}_methylation3", "motif_profile.csv")
        depth_file = os.path.join(work_dir, prefix, f"{prefix}_methylation3", "mean_depth.csv")
        fai = os.path.join(work_dir, prefix, f"{prefix}.hifiasm.p_ctg.rename.fa.fai")
        all_host_file = f"{all_dir}/{prefix}/all_host_ctgs.tsv" 
        ## skip if all_host_file does not exist
        if not os.path.exists(all_host_file):
            print(f"Skipping {prefix} as all_host_file does not exist.")
            continue
        if not os.path.exists(fai):
            continue
        best_ctgs, good_depth = get_best_ctg(depth_file, fai)
        best_ctgs = read_all_host(all_host_file, good_depth)
        ## read length of each ctg in a dict

        ## skip if not exist
        if not os.path.exists(drep_clu_file) or not os.path.exists(profile_file):
            print(f"Skipping {prefix} due to missing files.")
            continue

        multiple_strain_ctg, ctg_phylum, clu_num = read_drep_cluster(drep_clu_file, best_ctgs)
        if clu_num < 2:
            continue
        ## read motif profile file
        profiles = pd.read_csv(profile_file, index_col=0)
        ## binary it by 0.5
        profiles_bin = (profiles >= 0.4).astype(int)
        profiles_bin = profiles_bin.T
        
        # ctg_list = list(multiple_strain_ctg.keys())
        ctg_list = list(ctg_phylum.keys())
        # print (profiles_bin)
        for i in range(len(ctg_list)):
            for j in range(i + 1, len(ctg_list)):
                ctg1 = ctg_list[i]
                ctg2 = ctg_list[j]
                ## skip if ctg1 not in profiles_bin or ctg2 not in profiles_bin
                if ctg1 not in profiles_bin.index or ctg2 not in profiles_bin.index:
                    continue
                if ctg1 == ctg2:
                    continue
                # print (profiles_bin.loc[ctg1])
                intersection = np.sum((profiles_bin.loc[ctg1] & profiles_bin.loc[ctg2]).astype(int))
                union = np.sum((profiles_bin.loc[ctg1] | profiles_bin.loc[ctg2]).astype(int))
                if union == 0:
                    continue
                jaccard_similarity = intersection / union if union > 0 else 0
                if ctg1 in multiple_strain_ctg:
                    cluster1 = multiple_strain_ctg[ctg1]
                else:
                    cluster1 = ctg1
                if ctg2 in multiple_strain_ctg:
                    cluster2 = multiple_strain_ctg[ctg2]
                else:
                    cluster2 = ctg2

                if cluster1 != cluster2:
                    # Calculate Jaccard similarity
                    jaccard_similarities.append([jaccard_similarity, "diff", prefix])
                else:
                    jaccard_similarities.append([jaccard_similarity, "same", prefix])
        # if len(jaccard_similarities) > 5:
        #     break
    print (f"Total {len(jaccard_similarities)} Jaccard similarities calculated.")
    ## plot the Jaccard similarities
    df_jaccard = pd.DataFrame(jaccard_similarities, columns=['jaccard_similarity', 'cluster_type', 'sample'])
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df_jaccard, hue='cluster_type', x='sample', y='jaccard_similarity', palette='Set2')
    plt.xticks(rotation=90)
    plt.title('Jaccard Similarity of Motifs in Multiple Strain Clusters')
    plt.ylabel('Jaccard Similarity')
    plt.savefig("../../tmp/results/jaccard_similarity_batch.pdf", dpi=300, bbox_inches='tight')
    plt.clf()   
    
def get_phylum(best_ctgs):
    import re
    bin_phylum = {}
    df = pd.read_csv(gtdb_file, sep="\t")
    for index, row in df.iterrows():
        if re.search("Unclassified", row["classification"]):
            taxon = "Unclassified"
        else:
            # print (row["classification"])
            taxon = row["classification"].split(";")[1].strip()
        bin_phylum[row["user_genome"]] = taxon
    ctg_bin = {}
    with open(bin_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            ctg, bin_id = line.strip().split("\t")
            ctg_bin[ctg] = bin_id
    
    ctg_phylum = {}
    for ctg in best_ctgs:
        bin_id = ctg_bin[ctg]
        if bin_id in bin_phylum:
            ctg_phylum[ctg] = bin_phylum[bin_id]
        else:
            ctg_phylum[ctg] = "Unclassified"
    return (ctg_phylum)

def out_best_ctgs(ref, best_ref, best_ctgs):
    ## use biopython to extract the best contigs from the reference fasta file
    from Bio import SeqIO
    from Bio.SeqIO.FastaIO import SimpleFastaParser
    with open(ref, "r") as f_in, open(best_ref, "w") as f_out:
        for record in SeqIO.parse(f_in, "fasta"):
            if record.id in best_ctgs:
                SeqIO.write(record, f_out, "fasta")
    print(f"Extracted {len(best_ctgs)} best contigs to {best_ref}")

def out_best_ctgs2(ref, best_ref, best_ctgs, best_ctg_dir):
    ## use biopython to extract the best contigs from the reference fasta file
    from Bio import SeqIO
    from Bio.SeqIO.FastaIO import SimpleFastaParser    
    with open(ref, "r") as f_in:
        for record in SeqIO.parse(f_in, "fasta"):
            if record.id in best_ctgs:
                with open(os.path.join(best_ctg_dir, str(record.id) + ".fasta"), "w") as f_out:
                    SeqIO.write(record, f_out, "fasta")
    print(f"Extracted {len(best_ctgs)} best contigs to {best_ref}")

def count_motifs(depth_file, best_ctgs, work_dir, prefix, environment):
    ## read depth file
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    data = []
    for index, row in depth_df.iterrows():
        if row['depth'] >= 10:
            good_depth[row['contig']] = row['depth']
    best_depth_ctg = []
    for ctg in best_ctgs:
        if ctg in good_depth:
            best_depth_ctg.append(ctg)
    print(f"Total {len(best_depth_ctg)} contigs with depth >= 10 found.")
    has_motif_ctg_num = 0
    motif_num_list = []
    for ctg in best_depth_ctg:
        motif_file = os.path.join(work_dir,"motifs", f"{ctg}.motifs.csv")
        # print (motif_file)
        if os.path.exists(motif_file):
            df_motif = pd.read_csv(motif_file)
            ## only keep themotifs with fraction >= 0.4, and nDetected >=100
            unique_motifs_num, unique_motifs = get_unique_motifs(df_motif)
            if unique_motifs_num > 0:
                has_motif_ctg_num += 1
            motif_num_list.append(unique_motifs_num)
            data.append([prefix, unique_motifs_num, environment,ctg])
        else:
            print(f"No motifs found for {ctg}.")
    # print (f"Total {has_motif_ctg_num} contigs with motifs found in the best contigs with depth >= 10.")
    # print (has_motif_ctg_num/ len(best_depth_ctg) * 100, "% of the best contigs with depth >= 10 have motifs.")
    # ## print mean and median of motif numbers
    # print(f"Mean motif number: {np.mean(motif_num_list)}")
    # print(f"Median motif number: {np.median(motif_num_list)}")
    # print(f"Max motif number: {np.max(motif_num_list)}")
    # print(f"Min motif number: {np.min(motif_num_list)}")
    # ## plot the distribution of motif numbers
    # plt.figure(figsize=(10, 6))
    # plt.hist(motif_num_list, bins=50, color='blue', alpha=0.7)
    # plt.title('Distribution of Motif Numbers in Best Contigs')
    # plt.xlabel('Number of Motifs')
    # plt.ylabel('Frequency')
    # plt.savefig(os.path.join("../../tmp/results", f"motif_num_distribution_{prefix}.png"), dpi=300, bbox_inches='tight')
    return data

def count_MT_num(anno_file):
    
    df = pd.read_csv(anno_file, sep="\t")
    # Extract contig name from 'Gene' column and add as a new column
    df['contig'] = df['Gene'].apply(lambda x: '_'.join(str(x).split('_')[:-1]))  
    ctg_MT_num = {}
    for index, row in df.iterrows():
        if row['contig'] not in ctg_MT_num:
            ctg_MT_num[row['contig']] = set()
        if row['Gene type'] == 'MT':
            ctg_MT_num[row['contig']].add(row['HMM'])
    return ctg_MT_num

def plot_MT_motif():
    ctg_MT_num = count_MT_num(anno_file)
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    for index, row in depth_df.iterrows():
        if row['depth'] >= 10:
            good_depth[row['contig']] = row['depth']
    best_depth_ctg = []
    for ctg in best_ctgs:
        if ctg in good_depth:
            best_depth_ctg.append(ctg)
    print(f"Total {len(best_depth_ctg)} contigs with depth >= 10 found.")
    has_motif_ctg_num = 0
    motif_num_list = []
    data = []
    for ctg in best_depth_ctg:
        motif_file = os.path.join(work_dir,"motifs", f"{ctg}.motifs.csv")
        # print (motif_file)
        if os.path.exists(motif_file):
            df_motif = pd.read_csv(motif_file)
            ## only keep themotifs with fraction >= 0.4, and nDetected >=100
            unique_motifs_num, unique_motifs = get_unique_motifs(df_motif)
            data.append({
                'contig': ctg,
                'motif_num': unique_motifs_num,
                'MT_num': len(ctg_MT_num[ctg]) if ctg in ctg_MT_num else 0
            })
        else:
            print(f"No motifs found for {ctg}.")
    df_data = pd.DataFrame(data)
    ## plot dot plot with seaborn
    plt.figure(figsize=(7, 6))
    sns.scatterplot(data=df_data, x='MT_num', y='motif_num', palette='viridis', s=100)
    plt.title('MTase Number vs Motif Number in Best Contigs')
    plt.xlabel('Number of MTases')
    plt.ylabel('Number of Motifs')
    plt.savefig(os.path.join("../../tmp/results", "MT_motif_num_distribution.png"), dpi=300, bbox_inches='tight')

def read_all_host(all_host_file, good_depth):
    """
    Read all host contigs from a file.
    """
    all_host_ctgs = {}
    with open(all_host_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            ctg, ctg, domain = line.strip().split("\t")
            if ctg not in good_depth:
                continue
            all_host_ctgs[ctg] = domain
    return all_host_ctgs

def read_mapping(map_sum):
    """
    map_read_count,raw_read_count,ratio
    985092,1043059,0.9444259624815087
    """
    if not os.path.exists(map_sum):
        return 0
    ## get the ratio
    f = open(map_sum, "r")
    for line in f:
        if line.startswith("#"):
            continue
        map_read_count, raw_read_count, ratio = line.strip().split(",")
        print(f"Map read count: {map_read_count}, Raw read count: {raw_read_count}, Ratio: {ratio}")
    f.close()
    return float(ratio)

def read_host(host_sum_file):
    line_num = sum(1 for line in open(host_sum_file) if line.strip())
    if line_num < 2:
        print(f"Host summary file {host_sum_file} is empty or has only header.")
        return 0
    # Read the host summary file
    host_sum = pd.read_csv(host_sum_file)
    linkage_num = 0
    for index, row in host_sum.iterrows():
        # if row["pvalue"] >= 0.05:
        #     continue
        # if row['final_score'] <= 0.6:
        #     continue
        if row["specificity"] >= 0.01:
            continue
        if row['final_score'] <= 0.5:
            continue
        linkage_num += 1
    return linkage_num

def read_orphan(orphan_file):
    regulate_motif_num = 0
    if not os.path.exists(orphan_file):
        return regulate_motif_num
    df = pd.read_csv(orphan_file)
    regulate_motif_set = set()
    for index, row in df.iterrows():
        if row["significant"] != True:
            continue
        motif_tag = f"{row['motif']}_{row['exact_pos']}"
        regulate_motif_set.add(motif_tag)
    regulate_motif_num = len(regulate_motif_set)
    return regulate_motif_num

def count_modified_base(work_dir, prefix, best_ctgs, length_dict, env, score_cutoff = 30):
    base_data = []
    for ctg in best_ctgs:
        length = length_dict[ctg]
        gff = f"{work_dir}/gffs/{ctg}.reprocess.gff"
        modified_num = 0
        modified_motif_num = 0
        if os.path.exists(gff):
            for line in open(gff, "r"):
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if int(fields[5]) < score_cutoff:
                    continue
                modified_num += 1
                if re.search(";motif=", fields[8]):
                    modified_motif_num += 1
        modified_ratio = modified_num / length
        modified_motif_ratio = modified_motif_num/ length
        motif_ratio =  modified_motif_num / modified_num if modified_num > 0 else 0
        base_data.append([prefix, ctg, length, modified_num, modified_motif_num, modified_ratio, modified_motif_ratio, motif_ratio, env])
    return base_data

def count_all_motif_num():
    all_data = []
    all_base_data = []
    genome_data = []
    all_dir = "/home/shuaiw/borg/paper/run2/"
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        if prefix == "ocean_1":
            continue
        print (f"Processing {prefix}...")
        work_dir = f"{all_dir}/{prefix}/{prefix}_methylation3"
        fai = f"{all_dir}/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
        map_sum = f"{all_dir}/{prefix}/{prefix}.align.count.csv"
        all_host_file = f"{all_dir}/{prefix}/all_host_ctgs.tsv"
        depth_file = os.path.join(work_dir, "mean_depth.csv")
        host_sum_file = os.path.join(work_dir, "host_summary.csv")
        orphan_file = os.path.join(work_dir, "regulatory_motif_enrichment.csv")

        if not os.path.exists(fai):
            print(f"Skipping {prefix} as fai file does not exist.")
            continue
        ## skip if all_host_file does not exist
        if not os.path.exists(all_host_file):
            print(f"Skipping {prefix} as all_host_file does not exist.")
            continue

        map_ratio = read_mapping(map_sum)
        linkage_num =  read_host(host_sum_file)
        regulate_motif_num = read_orphan(orphan_file)
        N50, genome_size = count_N50_size(fai)
        print(f"{prefix}: N50 size: {N50}, Genome size: {genome_size}")

        good_depth, length_dict = get_best_ctg(depth_file, fai)
        best_ctgs = read_all_host(all_host_file, good_depth)
        print ("best ctgs num:", len(best_ctgs))

        genome_data.append([prefix, N50, genome_size, sample_env_dict[prefix], map_ratio, linkage_num, \
                            regulate_motif_num, len(best_ctgs)])

        # print (f"Total {len(best_ctgs)} best contigs with depth >= 10 found.")
        all_data += count_motifs(depth_file, best_ctgs, work_dir, prefix, sample_env_dict[prefix])
        all_base_data += count_modified_base(work_dir, prefix, best_ctgs, length_dict, sample_env_dict[prefix])
        # break
    print ("start plot...")
    df_all_data = pd.DataFrame(all_data, columns=['sample', 'motif_num', 'environment', 'contig'])
    df_genome_data = pd.DataFrame(genome_data, columns=['sample', 'N50', 'genome_size', 'environment', 'map_ratio', 'linkage_num', 'regulate_motif_num','best_ctg_num'])
    df_all_base_data = pd.DataFrame(all_base_data, columns=['sample', 'ctg', 'length', 'modified_num', 'modified_motif_num', 'modified_ratio', 'modified_motif_ratio', 'motif_ratio', 'environment'])
    plot_motif(df_all_data)
    plot_genome(df_genome_data)
    plot_meta(df_genome_data)
    plot_base(df_all_base_data)
    ## save genome data
    df_genome_data.to_csv("../../tmp/results2/genome_data_all_samples.csv", index = False)
    df_all_base_data.to_csv("../../tmp/results2/base_count_all_samples.csv", index=False)
    df_all_data.to_csv("../../tmp/results2/motif_num_all_samples.csv", index=False)
    
def get_stastics():
    df_genome_data = pd.read_csv("../../tmp/results2/genome_data_all_samples.csv")
    df_all_base_data = pd.read_csv("../../tmp/results2/base_count_all_samples.csv")
    df_all_data = pd.read_csv("../../tmp/results2/motif_num_all_samples.csv")
    ## exclude 96plex environment for df_genome_data, df_all_base_data, df_all_data
    df_genome_data = df_genome_data[df_genome_data['environment'] != '96plex']
    df_all_base_data = df_all_base_data[df_all_base_data['environment'] != '96plex']
    df_all_data = df_all_data[df_all_data['environment'] != '96plex']

    ## count sample number for each environment in df_genome_data
    env_groups = df_genome_data.groupby('environment')
    print ("Sample number for each environment:")
    for env, group in env_groups:
        total_samples_env = group.shape[0]
        print (f"Environment: {env}, Total samples: {total_samples_env}")
    print ("*************************")

    print ("Genome data summary:")
    ## count how many samples proportion have map_ratio > 0.5 and samples proportion have map_ratio > 0.7
    total_samples = df_genome_data.shape[0]
    map_ratio_50 = df_genome_data[df_genome_data['map_ratio'] > 0.5].shape[0]
    map_ratio_70 = df_genome_data[df_genome_data['map_ratio'] > 0.7].shape[0]
    print (f"Total samples: {total_samples}, Samples with map_ratio > 0.5: {map_ratio_50} ({map_ratio_50/total_samples*100:.2f}%), Samples with map_ratio > 0.7: {map_ratio_70} ({map_ratio_70/total_samples*100:.2f}%)")

    ## get proportion with motif_num > 0 in df_all_data
    total_ctgs = df_all_data.shape[0]
    ctgs_with_motif = df_all_data[df_all_data['motif_num'] > 0].shape[0]
    print (f"Total contigs: {total_ctgs}, Contigs with motif_num > 0: {ctgs_with_motif} ({ctgs_with_motif/total_ctgs*100:.2f}%)")
    ## count the average motif_num for all contigs
    avg_motif_num = df_all_data['motif_num'].mean()
    print (f"Average motif_num for all contigs: {avg_motif_num:.2f}") 
    ## count the proportion with motif_num > 0 for each environment
    env_groups = df_all_data.groupby('environment')
    for env, group in env_groups:
        total_ctgs_env = group.shape[0]
        ctgs_with_motif_env = group[group['motif_num'] > 0].shape[0]
        ## also print the average motif_num for each environment
        avg_motif_num_env = group['motif_num'].mean()
        print (f"Environment: {env}, Total contigs: {total_ctgs_env}, Contigs with motif_num > 0: {ctgs_with_motif_env}({ctgs_with_motif_env/total_ctgs_env*100:.2f}%), Average motif_num: {avg_motif_num_env:.2f}")
    print ("*************************")
    ## count average regulate_motif_num for each environment
    env_groups = df_genome_data.groupby('environment')
    for env, group in env_groups:
        total_samples_env = group.shape[0]
        avg_regulate_motif_num_env = group['regulate_motif_num'].mean()
        print (f"Environment: {env}, Total samples: {total_samples_env}, Average regulate_motif_num: {avg_regulate_motif_num_env:.2f}")
    ## count the proportion of infant sample has regulate_motif_num > 0
    infant_group = df_genome_data[df_genome_data['environment'] == 'infant']
    total_infant_samples = infant_group.shape[0]
    infant_samples_with_regulate_motif = infant_group[infant_group['regulate_motif_num'] > 0].shape[0]
    print (f"Infant samples: {total_infant_samples}, Infant samples with regulate_motif_num > 0: {infant_samples_with_regulate_motif} ({infant_samples_with_regulate_motif/total_infant_samples*100:.2f}%)")
    print ("*************************")
    ### count the average  modified_ratio for each environment
    env_groups = df_all_base_data.groupby('environment')
    for env, group in env_groups:
        total_ctgs_env = group.shape[0]
        avg_modified_ratio_env = group['modified_ratio'].mean()
        avg_modified_motif_ratio_env = group['modified_motif_ratio'].mean()
        avg_motif_ratio_env = group['motif_ratio'].mean()
        print (f"Environment: {env}, Total contigs: {total_ctgs_env}, Average modified_ratio: {avg_modified_ratio_env:.4f}, Average modified_motif_ratio: {avg_modified_motif_ratio_env:.4f}, Average motif_ratio: {avg_motif_ratio_env:.4f}")
    print ("*************************")

def plot_motif(df_all_data):
    ### sort by sample
    df_all_data = df_all_data.sort_values(by='sample')
    ## plot boxplot where sample is on x-axis and motif_num is on y-axis
    plt.figure(figsize=(14, 6))
    sns.boxplot(data=df_all_data, x='sample', y='motif_num', hue='environment', palette='Set2')
    plt.xticks(rotation=90)
    plt.title('Distribution of Motif Numbers Across Samples')
    plt.xlabel('Sample')
    plt.ylabel('Number of Motifs')
    plt.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.savefig("../../tmp/results2/motif_num_distribution_all_samples.png", dpi=300, bbox_inches='tight')

def plot_base(df_all_base_data):
    ### sort by sample
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 12))
    df_all_base_data = df_all_base_data.sort_values(by='sample')
    ## plot boxplot where sample is on x-axis and modified_ratio is on y-axis
    sns.boxplot(data=df_all_base_data, x='sample', y='modified_ratio', hue='environment', palette='Set2', ax=ax1)
    ax1.set_title('Distribution of Modified Ratios Across Samples')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Modified Ratio')
    ax1.tick_params(axis='x', rotation=90)
    ax1.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

    ## plot boxplot where sample is on x-axis and modified_motif_ratio is on y-axis
    sns.boxplot(data=df_all_base_data, x='sample', y='modified_motif_ratio', hue='environment', palette='Set2', ax=ax2)
    ax2.set_title('Distribution of Modified Motif Ratios Across Samples')
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Modified Motif Ratio')
    ax2.tick_params(axis='x', rotation=90)
    ax2.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

    ## plot boxplot where sample is on x-axis and motif_ratio is on y-axis
    sns.boxplot(data=df_all_base_data, x='sample', y='motif_ratio', hue='environment', palette='Set2', ax=ax3)
    ax3.set_title('Distribution of Motif Ratios Across Samples')
    ax3.set_xlabel('Sample')
    ax3.set_ylabel('Motif Ratio')
    ax3.tick_params(axis='x', rotation=90)
    ax3.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig("../../tmp/results2/base_count_all_samples.png", dpi=300, bbox_inches='tight')
    ### save the dataframe
    
def plot_meta(df_genome_data):
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 12))
    df_genome_data = df_genome_data.sort_values(by='sample')
    df_genome_data.set_index('sample', inplace=True)

    sns.barplot(x=df_genome_data.index, y='map_ratio', hue='environment', data=df_genome_data.reset_index(), ax=ax1, palette='Set2')
    ax1.set_title('Mapping Ratio Across Samples')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Mapping Ratio ')
    ax1.tick_params(axis='x', rotation=90)

    sns.barplot(x=df_genome_data.index, y='linkage_num', hue='environment', data=df_genome_data.reset_index(), ax=ax2, palette='Set2')
    ax2.set_title('Linkage Number Across Samples')
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Linkage Number ')
    ax2.tick_params(axis='x', rotation=90)

    sns.barplot(x=df_genome_data.index, y='regulate_motif_num', hue='environment', data=df_genome_data.reset_index(), ax=ax3, palette='Set2')
    ax3.set_title('Regulatory Motif Number Across Samples')
    ax3.set_xlabel('Sample')
    ax3.set_ylabel('Regulatory Motif Number ')
    ax3.tick_params(axis='x', rotation=90)

    ax1.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    # Remove legend for ax2 and ax3
    ax2.legend_.remove()
    ax3.legend_.remove()
    plt.tight_layout()
    plt.savefig("../../tmp/results2/meta_all_samples.png", dpi=300, bbox_inches='tight')  

def plot_genome(df_genome_data):
    ### plot bar plot for genome size and N50 in separate subplots, using log scale for y-axis
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20, 12))
    df_genome_data = df_genome_data.sort_values(by='sample')
    df_genome_data.set_index('sample', inplace=True)
    
    # N50 subplot
    sns.barplot(x=df_genome_data.index, y='N50', hue='environment', data=df_genome_data.reset_index(), ax=ax1, palette='Set2')
    ax1.set_yscale('log')
    ax1.set_title('N50 Across Samples')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('N50 Size (log scale)')
    ax1.tick_params(axis='x', rotation=90)
    ax1.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Genome size subplot
    sns.barplot(x=df_genome_data.index, y='genome_size', hue='environment', data=df_genome_data.reset_index(), ax=ax2, palette='Set2')
    ax2.set_yscale('log')
    ax2.set_title('Genome Size Across Samples')
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Genome Size (log scale)')
    ax2.tick_params(axis='x', rotation=90)
    ax2.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')

    sns.barplot(x=df_genome_data.index, y='best_ctg_num', hue='environment', data=df_genome_data.reset_index(), ax=ax3, palette='Set2')
    ax3.set_title('Best Contig Number Across Samples')
    ax3.set_xlabel('Sample')
    ax3.set_ylabel('Best Contig Number')
    ax3.tick_params(axis='x', rotation=90)
    ax3.legend(title='Environment', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig("../../tmp/results2/genome_size_N50_all_samples.png", dpi=300, bbox_inches='tight')

def count_N50_size(fai):
    """
    Count N50 size from a fasta index file.
    """
    total_length = 0
    contig_lengths = []
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            contig_lengths.append(length)
            total_length += length
    contig_lengths.sort(reverse=True)
    
    n50_size = 0
    half_length = total_length / 2
    current_length = 0
    for length in contig_lengths:
        current_length += length
        if current_length >= half_length:
            n50_size = length
            break
    return n50_size, total_length

def read_metadata(meta_file):
    sample_env_dict = {}
    with open(meta_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            sample = parts[1]
            env = parts[3]
            sample_env_dict[sample] = env
    return sample_env_dict

if __name__ == "__main__":
    meta_file = "/home/shuaiw/Methy/assembly_pipe/prefix_table.tab"
    sample_env_dict = read_metadata(meta_file)
    count_all_motif_num()
    # get_stastics()
    # jaccard()
    # jaccard_batch()

    # fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
    # ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
    # best_ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.circular.contigs.fa"
    # best_ctg_dir = "/home/shuaiw/methylation/data/borg/contigs/circular/"
    # profile_file = "/home/shuaiw/borg/bench/soil/run1/motif_profile.csv"
    # anno_file = "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_RM.rm.genes.tsv"
    # bin_file = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.bin.tab"
    # gtdb_file = "/home/shuaiw/borg/contigs/GTDB/gtdbtk.all.summary.tsv"

    # work_dir = "/home/shuaiw/borg/bench/soil/run1/"
    # fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
    # prefix = "soil"

    # work_dir = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2"
    # fai = "/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa.fai"
    # prefix = "cow"

    # work_dir = "/home/shuaiw/borg/paper/infant/NANO_2_INF1340011_4PB"
    # fai = "/home/shuaiw/borg/allison/NANO_2_INF1340011_4PB_HR_HIFIASM_META_scaffold_min1000.fa.fai"
    # prefix = "infant"
    # prefix = "SRR23446540"
    # work_dir = f"/home/shuaiw/borg/paper/run/{prefix}/{prefix}_methylation"
    # fai = f"/home/shuaiw/borg/paper/run/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
    # all_host_file = f"/home/shuaiw/borg/paper/run/{prefix}/all_host_ctgs.tsv"


    # jaccard()
    # plot_MT_motif()
    # out_best_ctgs2(ref, best_ref, best_ctgs, best_ctg_dir)  ## split ctgs


    # out_best_ctgs(ref, best_ref, best_ctgs)
    # plot_heatmap(profile_file)
    # plot_strain_heatmap(profile_file)