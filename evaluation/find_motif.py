import numpy as np
import pandas as pd
from collections import defaultdict
import subprocess
import logging
import os
import random


from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET

METHYLATION_TYPES = {"6mA": "a", "5mC": "m", "4mC": "21839", "5hmC": "h"}
METHYLATION_TYPES_REV = {"a": "6mA", "m": "5mC", "21839": "4mC", "h": "5hmC"}
REF = {}
WINDOW_SIZE = 12
MIN_EVALUE = 0.1
MOTIF_FREQ_CUTOFF = (
    0.8  # If a nucleotide isn't 80% of sites in a motif, converted to an N
)

def check_existing_motif(attributes_dict):
    """
    check if motif maker has detected the motif
    """
    if 'motif' not in attributes_dict:
        motif = 'None'
    else:
        motif = attributes_dict['motif']
    return motif

def run_streme(kmer_file, streme_path="streme"):
    """Run STREME on a set of files
    Args:
        kmer_file: The FASTA file of methylated kmers, returned by write_to_fasta
        streme_path: Path to the STREME program.
    Returns:
        output_dir: path to output directory of STREME
    """

    ## Run STREME
    output_dir = kmer_file.split(".fasta")[0] + "_streme"
    ## if the output directory exists, remove it
    if os.path.exists(output_dir):
        os.system(f"rm -rf {output_dir}")
    print ("output_dir", output_dir)

    cmd = f"streme .  -p {kmer_file} -oc {output_dir} --dna -n /home/shuaiw/borg/bench/test/control.fasta"
    # cmd = f"meme {kmer_file} -oc {output_dir} -minw 2 -maxw 10 -dna -revcomp -p 10"
    

    # cmd = cmd.format(
    #     streme=streme_path,
    #     output=output_dir,
    #     input_file=kmer_file,
    # )
    logging.info("Running STREME: %s", cmd)
    print (cmd)

    subprocess.check_output(
        cmd, shell=True, stderr=subprocess.DEVNULL, universal_newlines=True
    )
    return output_dir

def read_streme(streme_output):
    """Read the output of STREME and assign methylated sites to any identified significant motifs.
    Args:
        modkit_table: The output of read_modkit()
        streme_output: The directory output of STREME
    Returns:
        modkit_table2: A modified version of modkit_table with the motif column added for motif assignments.
    """

    tree = ET.parse(streme_output + "/streme.xml")
    root = tree.getroot()

    motifs = []
    motif_new_to_original = {}  # store original motif names
    motif_sites = {}

    for motif in root[1]:

        evalue = float(motif.get("test_evalue"))
        print ("evalue", evalue)
        if evalue < MIN_EVALUE:
            
            ### Evalue cutoff
            motif_old = motif.get("id").split("-")[1]
            motif_new = ""

            ### Convert to N's
            i = 0
            for pos in motif:
                freqs = [float(x) for x in pos.attrib.values()]
                if motif_old[i] in ["A", "C", "G", "T"]:  # If it's a single base
                    if max(freqs) <= MOTIF_FREQ_CUTOFF:  # If it's not >80%, make it N
                        motif_new += "N"
                    else:
                        motif_new += motif_old[i]
                elif motif_old[i] in ["W", "S", "M", "K", "R", "Y"]:  # If it's a double
                    top_two = sum(sorted(freqs)[-2:])
                    if top_two <= MOTIF_FREQ_CUTOFF:  # If it's not >80%, make it N
                        motif_new += "N"
                    else:
                        motif_new += motif_old[i]
                elif motif_old[i] in ["B", "D", "H", "V"]:  # All triplies become N
                    motif_new += "N"
                else:
                    motif_new += motif_old[i]
                i += 1

            ### trim N's
            motif_new = motif_new.strip("N")
            print(
                "Found and trimmed motif: "
                + motif.get("id")
                + "\t"
                + motif_new
                + "\t"
                + motif.get("total_sites")
            )
            motifs.append(motif_new)
            motif_new_to_original[motif_new] = motif.get("id")

            ## Find motif occurrences in reference
            motif_len = len(motif_new)

            for r, contig in REF.items():
                for site in nt_search(str(contig), motif_new)[1:]:
                    for i in range(site, site + motif_len):
                        motif_sites[r + ":" + str(i) + "+"] = motif_new

                for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
                    1:
                ]:
                    for i in range(site, site + motif_len):
                        motif_sites[r + ":" + str(i) + "-"] = motif_new
    print ("motif number", len(motifs), motifs, len(motif_sites))
    # Print the first 5 elements of motif_sites
    # for i, (key, value) in enumerate(motif_sites.items()):
    #     if i >= 100:
    #         break
    #     print(f"{key}: {value}")
    return motifs, motif_sites
# Alternatively, you can parse the GFF file manually without using gffutils
def parse_gff(file_path):
    motif_info_dict = defaultdict(dict)
    motif_dict = {}
    context_f = open(context_fasta, 'w')
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            attributes_dict = {key: value for key, value in (item.split('=') for item in attributes.split(';'))}

            print (f">{seqid}_{start}\t{check_existing_motif(attributes_dict)}\n{attributes_dict['context']}", file = context_f)

            if 'motif' not in attributes_dict:
                ### check reason for this
                continue

            # print (attributes_dict)
            if attributes_dict['motif'] not in motif_info_dict[seqid]:
                motif_info_dict[seqid][attributes_dict['motif']] = []
            
            # motif_info_dict[seqid][attributes_dict['motif']].append(float(attributes_dict['IPDRatio']))
            # motif_info_dict[seqid][attributes_dict['motif']].append(int(score))
            # motif_info_dict[seqid][attributes_dict['motif']].append(float(attributes_dict['frac']))

            if int(score) >= 30:
                motif_info_dict[seqid][attributes_dict['motif']].append(1)
            else:
                motif_info_dict[seqid][attributes_dict['motif']].append(0)
            
            if attributes_dict['motif'] not in motif_dict:
                motif_dict[attributes_dict['motif']] = 0
            motif_dict[attributes_dict['motif']] += 1
    print ("motif number", len(motif_dict), 'seq number', len(motif_info_dict))
    context_f.close()
    return motif_info_dict, motif_dict

def add_benchmark(motif_info_dict, motif_dict):
    new_motif_info_dict = defaultdict(dict)
    for seqid in motif_info_dict:
        new_motif_info_dict[seqid] = motif_info_dict[seqid]
        for motif in motif_dict:
            new_motif_info_dict[seqid][motif] = motif_info_dict[seqid][motif]
            if motif in motif_info_dict[seqid]:
                motif_score_num = len(motif_info_dict[seqid][motif])
                half = int(motif_score_num/2)
                for i in range(2):
                    new_seq_id = seqid + "_" + str(i)
                    new_score_list = motif_info_dict[seqid][motif][half*i:half*(i+1)]
                    # print ("new_seq_id", new_seq_id, "motif", motif, "score", new_score_list[:10])
                    new_motif_info_dict[new_seq_id][motif] = new_score_list
    # print (new_motif_info_dict.keys())
    return new_motif_info_dict

#get motif matrix
def get_motif_matrix(motif_info_dict, motif_dict):
    data = []
    for motif in motif_dict:
        for seqid in motif_info_dict:
            # print(f"Sequence ID: {seqid}")
            if motif in motif_info_dict[seqid]:
                # print (motif_info_dict[seqid][motif])
                average_score = np.mean(motif_info_dict[seqid][motif])
                # print ("motif", motif, "seqid", seqid, "score", motif_info_dict[seqid][motif][:10])
                # average_score = np.median(motif_info_dict[seqid][motif])
            else:
                average_score = 0
            data.append([seqid, motif, average_score])
    # for da in data:
    #     print(da)   
    df = pd.DataFrame(data, columns = ['seqid', 'motif', 'score'])
    ## output to csv
    df.to_csv(matrix_file, index = False)

def heatmap(matrix_file):
    ## set figure size
    
    ## given dataframe, plot the heatmap
    df = pd.read_csv(matrix_file)
    df = df.pivot(index="seqid", columns="motif", values="score")
    # df = df.pivot("seqid", "motif", "score")
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.figure(figsize=(30, 15))
    # sns.heatmap(df)
    # Plot the heatmap with hierarchical clustering
    sns.clustermap(df, method='average', metric='euclidean', cmap='viridis', figsize=(30, 15))

    
    plt.savefig(figure_name)

def cluster_contigs(matrix_file):
    ## given dataframe, plot the heatmap
    df = pd.read_csv(matrix_file)
    df = df.pivot(index="seqid", columns="motif", values="score")
    ## get the sequence id list
    seqid_list = df.index.tolist()
    ## original contig number list
    raw_seqid_list = seqid_list
    for seqid in seqid_list:
        raw_seqid_list.append("_".join(seqid.split("_")[:-2]))
    ## covert df to matrix, and then to array
    matrix = df.to_numpy()
    print (matrix.shape)
    ## reduce dimention using t-SNE
    from sklearn.manifold import TSNE
    X_embedded = TSNE(n_components=2).fit_transform(matrix)
    print (X_embedded)
    print (X_embedded.shape)
    ## plot the t-SNE result
    import matplotlib.pyplot as plt
    plt.scatter(X_embedded[:, 0], X_embedded[:, 1])
    ## save the figure in pdf
    plt.savefig("break_contigs/break_contigs.tSNE.pdf")
    ## cluster the contigs using DBSCAN based on the t-SNE result
    from sklearn.cluster import DBSCAN
    clustering = DBSCAN(eps=0.5, min_samples=2).fit(X_embedded)
    print (clustering.labels_)
    ## plot the cluster result
    plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c = clustering.labels_, shape = raw_seqid_list)
    ## save the figure in pdf
    plt.savefig("break_contigs/break_contigs.cluster.pdf")

def separat_motif(raw_gff):
    motif_info_dict = defaultdict(dict)
    motif_dict = {}
    context_fasta = raw_gff[:-4] + ".context.fasta"
    context_f = open(context_fasta, 'w')
    i = 1
    context_dict = {}
    with open(raw_gff, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            attributes_dict = {key: value for key, value in (item.split('=') for item in attributes.split(';'))}
            if float(attributes_dict['IPDRatio']) > 3:
                continue
            context_dict[attributes_dict['context']] = float(attributes_dict['IPDRatio'])

    ### sort the context_dict by IPDRatio
    context_dict = dict(sorted(context_dict.items(), key=lambda item: item[1], reverse=True))
    for context in context_dict:
        if i > 0:
            print (context, context_dict[context])
            print (f">context_{i}\n{context}", file = context_f)
            # print (f">{seqid}_{start}\t{check_existing_motif(attributes_dict)}\n{attributes_dict['context']}", file = context_f)
        if i == 10000:
            break
        i += 1
    context_f.close()
    run_streme(context_fasta)

def merge_gff(gff_folder, gff_file):
    os.system(f"cat {gff_folder}/*.reprocess.gff > {gff_file}")

def load_motifs(folder):
    motif_dict = {}
    motif_info_dict = defaultdict(dict)
    for file in os.listdir(folder):
        if file.endswith(".motif.csv"):
            ## get base name of the file
            seq_id = file.split(".motif.csv")[0]
            ### read the motif file using pandas
            df = pd.read_csv(folder + file)
            for index, row in df.iterrows():
                motif = row['motifString']
                motif_info_dict[seq_id][motif] = float(row['fraction'])
                motif_dict[motif] = 1
    print ("motif number", len(motif_dict), 'seq number', len(motif_info_dict))
    return motif_info_dict, motif_dict

def get_reverse_complement_seq(sequence):
    sequence = sequence[::-1]
    trantab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
    string = sequence.translate(trantab)
    return string

def load_ipd_ratio(csv, ref_seq):
    df = pd.read_csv(csv, sep = ",")
    for_dict = {}
    rev_dict = {}
    for_dict_all = {}
    rev_dict_all = {}
    for index, row in df.iterrows():
        if row['strand'] == 1:
            for_dict_all[row['tpl']] = [row['ipd_ratio'], row['pvalue']]
        else:
            rev_dict_all[row['tpl']] = [row['ipd_ratio'], row['pvalue']]
    ## filter the ipd_ratio with pvalue < 0.01
    df = df[df['pvalue'] < 0.01]
    for index, row in df.iterrows():
        if row['strand'] == 1:
            for_dict[row['tpl']] = row['pvalue']
        else:
            rev_dict[row['tpl']] = row['pvalue']
    
    context_f = open(context_fasta, 'w')
    for site in for_dict:
        start = site - flank_len
        end = site + flank_len
        if start >= 0 and end <= len(ref_seq):
            seq = ref_seq[start:end]
            # print (f">for_{site}\n{seq}")
            print (f">for_{site}\n{seq}", file = context_f)

    for site in rev_dict:
        start = site - flank_len
        end = site + flank_len
        if start >= 0 and end <= len(ref_seq):
            seq = ref_seq[start:end]
            # print (f">rev_{site}\n{get_reverse_complement_seq(seq)}")
            print (f">rev_{site}\n{get_reverse_complement_seq(seq)}", file = context_f)
    context_f.close()

    return for_dict, rev_dict, for_dict_all, rev_dict_all

def generate_control(control_fasta, ref_seq, for_dict, rev_dict, max_seq_num, flank_len):
    ### randomly select max_seq_num sites from ref_seq and the site should not in for_dict and rev_dict
    control_f = open(control_fasta, 'w')
    seq_len = len(ref_seq)
    sites = []
    for i in range(max_seq_num):
        site = random.randint(0, seq_len)
        while site in for_dict or site in rev_dict:
            site = random.randint(0, seq_len)
        sites.append(site)
    ### print the flank sequence of the sites from ref
    for site in sites:
        start = site - flank_len
        end = site + flank_len
        if start >= 0 and end <= seq_len:
            seq = ref_seq[start:end]
            # print (site, ref_seq[start:end])
            ## 50 % percent should be reverse complement
            if random.random() > 0.5:
                seq = get_reverse_complement_seq(seq)
            print (f">control_{site}\n{seq}", file = control_f)
    control_f.close()
    
def read_ref(ref):
    seq_dict = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        return str(record.seq), record.id

def motif_stastics(ref_id, motifs, motif_sites, for_dict, rev_dict):
    for i in for_dict:
        tag = ref_id + ":" + str(i) + "+"
        # print (tag)
        if tag in motif_sites:
            print (tag, motif_sites[tag], for_dict[i])

### SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L borg
### SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_494_C  Methanoperedens_43_0_circular
#### SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_267_C Methanoperedens_43_25_c

if __name__ == "__main__":

    # raw_gff = "/home/shuaiw/methylation/data/borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_10_L.gff"
    # # raw_gff = "/home/shuaiw/methylation/data/borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_11_C.gff"
    # ref="/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa"
    # ipd_ratio ="/home/shuaiw/methylation/data/borg/b_contigs/test2/ipd_ratio/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_0_219069.ipd3.csv"
    # streme_output = "/home/shuaiw/borg/streme_output"
    # context_fasta = "/home/shuaiw/methylation/data/borg/b_contigs/test2/" + "context.fasta"
    # control_fasta = "/home/shuaiw/methylation/data/borg/b_contigs/test2/" + "control.fasta"


    # raw_gff = "/home/shuaiw/borg/bench/test/E_coli_H10407_1.gff"
    raw_gff = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/gffs/E_coli_H10407_1.gff"
    # raw_gff = "/home/shuaiw/borg/bench/zymo_new_ref/gffs/E_coli_H10407_1.gff"
    ref = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/contigs/E_coli_H10407_1.fa"
    # ipd_ratio = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/ipd_ratio/E_coli_H10407_1.ipd3.csv"
    # streme_output = "/home/shuaiw/borg/bench/test/streme_output"
    # context_fasta = "/home/shuaiw/borg/bench/test/context.fasta"
    # control_fasta = "/home/shuaiw/borg/bench/test/control.fasta"


    max_seq_num = 1000
    flank_len = 25
    ref_seq, ref_id = read_ref(ref)
    # for_dict, rev_dict, for_dict_all, rev_dict_all = load_ipd_ratio(ipd_ratio, ref_seq)
    # generate_control(control_fasta, ref_seq, for_dict, rev_dict, max_seq_num, flank_len)
    # read_streme(streme_output)

    separat_motif(raw_gff)

    # motifs, motif_sites = read_streme (streme_output)
    # motif_stastics(ref_id, motifs, motif_sites, for_dict_all, rev_dict_all)

    # streme -p /home/shuaiw/methylation/data/borg/b_contigs/test2/context.fasta -n /home/shuaiw/methylation/data/borg/b_contigs/test2/control.fasta --oc ~/borg/streme_output