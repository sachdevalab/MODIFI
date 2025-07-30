import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics.pairwise import cosine_similarity

from Bio.Seq import Seq


def get_best_ctg(depth_file, fai, min_len = 1000000):
    """
    Get the best contig based on length from a fasta file.
    """
    depth_df = pd.read_csv(depth_file)
    good_depth = {}
    for index, row in depth_df.iterrows():
        if row['depth'] >= 10:
            good_depth[row['contig']] = row['depth']
    best_ctgs = []
    with open(fai, "r") as f:
        for line in f:
            ctg, length, _, _, _ = line.strip().split("\t")
            length = int(length)
            if ctg[-1] == "C" and length >= min_len and ctg in good_depth:
                best_ctgs.append(ctg)
    print (f"Total {len(best_ctgs)} contigs with length >= {min_len} found.")
    return best_ctgs, good_depth

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

def read_drep_cluster(drep_clu_file):
    
    drep_clu_dict = {}
    df = pd.read_csv(drep_clu_file)
    for index, row in df.iterrows():
        drep_clu_dict[row['genome'][:-6]] = row['secondary_cluster']
    ## count the number of contigs in each cluster
    clu_count = {}
    for ctg, clu in drep_clu_dict.items():
        if clu not in clu_count:
            clu_count[clu] = 0
        clu_count[clu] += 1
    multiple_strain_ctg = {}
    ctg_phylum = {}
    for ctg, clu in drep_clu_dict.items():
        if clu_count[clu] > 1:
            multiple_strain_ctg[ctg] = clu
            ctg_phylum[ctg] = clu
    print ("multiple_strain_ctg", multiple_strain_ctg)
    return multiple_strain_ctg, ctg_phylum

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

def get_unique_motifs(df_motif):
    df_motif = df_motif[(df_motif['fraction'] >= 0.4) & (df_motif['nDetected'] >= 100)]
    ## rm redundant motifs which are reverse complement 
    unique_motifs = []
    for index, row in df_motif.iterrows():
        if row['motifString'] not in unique_motifs and  str(Seq(row['motifString']).reverse_complement()) not in unique_motifs:
            unique_motifs.append(row['motifString'])
    return len(unique_motifs)

def count_motifs(depth_file, best_ctgs, work_dir, prefix):
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
            unique_motifs_num = get_unique_motifs(df_motif)
            if unique_motifs_num > 0:
                has_motif_ctg_num += 1
            motif_num_list.append(unique_motifs_num)
            data.append([prefix,unique_motifs_num ])
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
            unique_motifs_num = get_unique_motifs(df_motif)
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

def count_all_motif_num():
    all_data = []
    all_dir = "/home/shuaiw/borg/paper/run2/"
    for my_dir in os.listdir(all_dir):
        prefix = my_dir
        print (f"Processing {prefix}...")
        work_dir = f"{all_dir}/{prefix}/{prefix}_methylation"
        fai = f"{all_dir}/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
        all_host_file = f"{all_dir}/{prefix}/all_host_ctgs.tsv"
        ## skip if all_host_file does not exist
        if not os.path.exists(all_host_file):
            print(f"Skipping {prefix} as all_host_file does not exist.")
            continue
        depth_file = os.path.join(work_dir, "mean_depth.csv")
        best_ctgs, good_depth = get_best_ctg(depth_file, fai)
        best_ctgs = read_all_host(all_host_file, good_depth)
        print (f"Total {len(best_ctgs)} best contigs with depth >= 10 found.")
        sample_data = count_motifs(depth_file, best_ctgs, work_dir, prefix)
        if not sample_data:
            print(f"No motifs found for {prefix}.")
            continue
        all_data += sample_data

    ## convert to df
    df_all_data = pd.DataFrame(all_data, columns=['sample', 'motif_num'])
    ## plot boxplot where sample is on x-axis and motif_num is on y-axis
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df_all_data, x='sample', y='motif_num', palette='Set2')
    plt.xticks(rotation=90)
    plt.title('Distribution of Motif Numbers Across Samples')
    plt.xlabel('Sample')
    plt.ylabel('Number of Motifs')
    plt.savefig("../../tmp/results/motif_num_distribution_all_samples.png", dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    count_all_motif_num()
    # jaccard()

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