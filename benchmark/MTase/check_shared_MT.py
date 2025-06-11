import pandas as pd
from collections import defaultdict
from Bio import pairwise2


def count_MT_num(anno_file):
    ctg_HMM_genes = defaultdict(dict)
    df = pd.read_csv(anno_file, sep="\t")
    # Extract contig name from 'Gene' column and add as a new column
    df['contig'] = df['Gene'].apply(lambda x: '_'.join(str(x).split('_')[:-1]))  
    ctg_MT_num = {}
    for index, row in df.iterrows():
        if row['contig'] not in ctg_MT_num:
            ctg_MT_num[row['contig']] = set()
        if row['Gene type'] == 'MT':
            ctg_MT_num[row['contig']].add(row['HMM'])
        if row['HMM'] not in ctg_HMM_genes[row['contig']]:
            ctg_HMM_genes[row['contig']][row['HMM']] = set()
        ctg_HMM_genes[row['contig']][row['HMM']].add(row['Gene'])
    return ctg_MT_num, ctg_HMM_genes

def read_protein(protein_file):
    ctg_protein = defaultdict(dict)
    ## read it using biopython
    from Bio import SeqIO
    for record in SeqIO.parse(protein_file, "fasta"):
        ctg_protein[record.id] = str(record.seq)
    return ctg_protein

def compare_ctgs_MTase(MGE_ctg, host_ctg, ctg_MT_num, ctg_HMM_genes, ctg_protein):
    shared_MT = ctg_MT_num[MGE_ctg].intersection(ctg_MT_num[host_ctg])
    # print ("shared MTase genes:", shared_MT)
    for mt in shared_MT:
        # print (mt, ctg_HMM_genes[MGE_ctg][mt], ctg_HMM_genes[host_ctg][mt])
        first_MGE_mt = list(ctg_HMM_genes[MGE_ctg][mt])[0]
        first_host_mt = list(ctg_HMM_genes[host_ctg][mt])[0]
        # print (ctg_protein[first_MGE_mt])
        # print (ctg_protein[first_host_mt])
        identity_list = []
        for first_MGE_mt in ctg_HMM_genes[MGE_ctg][mt]:
            for first_host_mt in ctg_HMM_genes[host_ctg][mt]:
                ## check the sequence similarity
                
                alignments = pairwise2.align.globalxx(ctg_protein[first_MGE_mt], ctg_protein[first_host_mt])
                for alignment in alignments:
                    identity = alignment[2] / alignment[4]
                    # print (mt, first_MGE_mt, first_host_mt)
                    # print(pairwise2.format_alignment(*alignment, full_sequences=True))
                    # print("Length of alignment:", alignment[4])  # Length of the alignment
                    # print("Score of alignment:", alignment[2])  # Score of the alignment
                    # print("Identity:", identity)  # Identity of the alignment
                    # print("\n")
                    identity_list.append(identity)
                    break
                # print ("-"*50)
                # print ("="*50)
        print (MGE_ctg, host_ctg, mt, "Max identity:", max(identity_list))

plasmid_list = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list"
rm_anno = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_RM.rm.genes.tsv"
protein_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_RM.rm.genes.faa"
ctg_MT_num, ctg_HMM_genes = count_MT_num(rm_anno)
ctg_protein = read_protein(protein_file)
MGE_ctg = "E_coli_H10407_2"
host_ctg = "E_coli_H10407_1"

## check the shared MTase genes
# print (ctg_MT_num[MGE_ctg])
# print (ctg_MT_num[host_ctg])

for line in open(plasmid_list):
    MGE_ctg = line.strip()
    if MGE_ctg == "seq_name":
        continue
    host_ctg = MGE_ctg[:-1] + "1"
    if host_ctg not in ctg_MT_num or MGE_ctg not in ctg_MT_num:
        continue
    compare_ctgs_MTase(MGE_ctg, host_ctg, ctg_MT_num, ctg_HMM_genes, ctg_protein)


