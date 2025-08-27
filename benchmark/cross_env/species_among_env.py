import pandas as pd
import re
from Bio import SeqIO
import sys
import os
from collections import defaultdict

def classify(classification):
    if re.search(r'd__Bacteria', classification) or re.search(r'Unclassified Bacteria', classification):
        return 'Bacteria'
    elif re.search(r'd__Archaea', classification) or re.search(r'Unclassified Archaea', classification):
        return 'Archaea'
    else:
        return classification

def get_bacteria_archea(gtdb_ar53, gtdb_bac120, depth_dict, min_dp = 10):
    classification_dict = defaultdict(list)
    ## check if the files exist
    if os.path.exists(gtdb_ar53):
        df_ar53 = pd.read_csv(gtdb_ar53, sep="\t")
        for index, row in df_ar53.iterrows():
            if row['user_genome'] not in depth_dict:
                continue
            if depth_dict[row['user_genome']] < min_dp:
                continue
            if row['classification'] == 'Unclassified':
                continue
            field = row['classification'].split(";")
            if len(field) < 7:
                continue
            species = field[6]
            classification_dict[species].append(row['user_genome'])

    else:
        print(f"Warning: {gtdb_ar53} does not exist. Skipping GTDB AR53 classification.")
    if os.path.exists(gtdb_bac120):
        df_bac120 = pd.read_csv(gtdb_bac120, sep="\t")
        for index, row in df_bac120.iterrows():
            if row['user_genome'] not in depth_dict:
                continue
            if depth_dict[row['user_genome']] < min_dp:
                continue
            if row['classification'] == 'Unclassified':
                continue
            field = row['classification'].split(";")
            if len(field) < 7:
                continue
            species = field[6]
            classification_dict[species].append(row['user_genome'])
    else:
        print(f"Warning: {gtdb_bac120} does not exist. Skipping GTDB BAC120 classification.")
    # print(classification_dict)
    return classification_dict


def read_depth(depth_file):
    depth_dict = {}
    df = pd.read_csv(depth_file, sep=",")
    for index, row in df.iterrows():
        depth_dict[row['contig']] = row['depth']
    return depth_dict

def each_sample(workdir, prefix):
    fai = f"{workdir}/{prefix}.hifiasm.p_ctg.rename.fa.fai"
    checkM_report = f"{workdir}/checkM2/quality_report.tsv"
    high_quality_file = f"{workdir}/all_host_ctgs.tsv"
    gtdb_ar53 = f"{workdir}/GTDB/gtdbtk.ar53.summary.tsv"
    gtdb_bac120 = f"{workdir}/GTDB/gtdbtk.bac120.summary.tsv"
    depth_file = f"{workdir}/{prefix}_methylation2/mean_depth.csv"
    depth_dict = read_depth(depth_file)
    classification_dict = get_bacteria_archea(gtdb_ar53, gtdb_bac120, depth_dict)
    print ("species number:", len(classification_dict))
    return classification_dict

def get_shared_species(classification_dict_1, classification_dict_2):
    shared_species = set(classification_dict_1.keys()).intersection(set(classification_dict_2.keys()))
    return shared_species

def get_motif(workdir, prefix, contig_set):
    motif_dict = defaultdict(list)
    for contig in contig_set:
        motif_file = f"{workdir}/{prefix}_methylation2/motifs/{contig}.motifs.csv"
        # print (motif_file)
        if os.path.exists(motif_file):
            df = pd.read_csv(motif_file, sep=",")
            for index, row in df.iterrows():
                motif_dict[contig].append(row['motifString'] + "_" + str(row["centerPos"]))
        else:
            print(f"Warning: {motif_file} does not exist.")
    return motif_dict

if __name__ == "__main__":
    # workdir = sys.argv[1]
    # prefix = sys.argv[2]
    
    prefix = "infant_3"
    workdir = f"/home/shuaiw/borg/paper/run2/{prefix}/"
    classification_dict_1 = each_sample(workdir, prefix)

    
    prefix2= "infant_4"
    workdir2 = f"/home/shuaiw/borg/paper/run2/{prefix2}/"
    classification_dict_2 = each_sample(workdir2, prefix2)

    shared_species = get_shared_species(classification_dict_1, classification_dict_2)
    print ("shared species:", shared_species, len(shared_species))

    ## for each shared species, check they the contig in each sample, and read the motif for each contig
    for species in shared_species:
        contigs_set1 = classification_dict_1[species]
        contigs_set2 = classification_dict_2[species]

        motif_dict_1 = get_motif(workdir, prefix, contigs_set1)
        print (species)
        print(prefix, motif_dict_1)
        motif_dict_2 = get_motif(workdir2, prefix2, contigs_set2)
        print(prefix2, motif_dict_2)
        # break

        print ("**********************\n")
