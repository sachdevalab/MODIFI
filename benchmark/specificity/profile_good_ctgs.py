import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE

from Bio.Seq import Seq


def get_best_ctg(min_len = 1000000):
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
    return best_ctgs

def plot_heatmap(profile_file, min_frac=0.4):
    profiles = pd.read_csv(profile_file, index_col=0)
    ## only keep the columns that are in the best contigs
    best_ctgs = get_best_ctg()
    profiles = profiles.loc[:, profiles.columns.isin(best_ctgs)]
    print ("original shape", profiles.shape)                                    
    # profiles = profiles.loc[(profiles > min_frac).any(axis=1)]
    # print ("filtered shape 1", profiles.shape)
    # Filter columns where any value is greater than 0.5
    # profiles = profiles.loc[:, (profiles > min_frac).any(axis=0)]

    # print ("filtered shape", profiles.shape)

    df = profiles.T

    sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(20, 15))
    plt.savefig("../../tmp/results/heatmap.png", dpi=300, bbox_inches='tight')
    plt.clf()

    drep_sim_dict = get_ctg_sim()
    ### for each contig in df, calculate the similarity with other contigs, and caculate the euclidean distance
    print (df)

def get_ctg_sim():
    drep_sim_dict = {}
    drep_sim = "/home/shuaiw/borg/contigs/dRep/dRep_out/data_tables/Mdb.csv"
    df = pd.read_csv(drep_sim)
    for index, row in df.iterrows():
        ctg1 = row['genome1'][:-5]
        ctg2 = row['genome2'][:-5]
        drep_sim_dict[ctg1 + "&" + ctg2] = row['similarity']
        drep_sim_dict[ctg2 + "&" + ctg1] = row['similarity']
    return drep_sim_dict





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

def count_motifs(depth_file, best_ctgs, work_dir):
    ## read depth file
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

        else:
            print(f"No motifs found for {ctg}.")
    print (f"Total {has_motif_ctg_num} contigs with motifs found in the best contigs with depth >= 10.")
    print (has_motif_ctg_num/ len(best_depth_ctg) * 100, "% of the best contigs with depth >= 10 have motifs.")
    ## print mean and median of motif numbers
    print(f"Mean motif number: {np.mean(motif_num_list)}")
    print(f"Median motif number: {np.median(motif_num_list)}")
    print(f"Max motif number: {np.max(motif_num_list)}")
    print(f"Min motif number: {np.min(motif_num_list)}")
    ## plot the distribution of motif numbers
    plt.figure(figsize=(10, 6))
    plt.hist(motif_num_list, bins=50, color='blue', alpha=0.7)
    plt.title('Distribution of Motif Numbers in Best Contigs')
    plt.xlabel('Number of Motifs')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join("../../tmp/results", "motif_num_distribution.png"), dpi=300, bbox_inches='tight')

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

fai = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa.fai"
ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
best_ref = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.circular.contigs.fa"
best_ctg_dir = "/home/shuaiw/methylation/data/borg/contigs/circular/"
profile_file = "/home/shuaiw/borg/bench/soil/run1/motif_profile.csv"
anno_file = "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_RM.rm.genes.tsv"

work_dir = "/home/shuaiw/borg/bench/soil/run1/"
depth_file = os.path.join(work_dir, "mean_depth.csv")

# best_ctgs = get_best_ctg()
# count_motifs(depth_file, best_ctgs, work_dir)
# plot_MT_motif()
# out_best_ctgs2(ref, best_ref, best_ctgs, best_ctg_dir)  ## split ctgs


# out_best_ctgs(ref, best_ref, best_ctgs)
plot_heatmap(profile_file)