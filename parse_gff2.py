import numpy as np
import pandas as pd
from collections import defaultdict
import subprocess
import logging
import os


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

    cmd = "{streme} .  -p {input_file} -oc {output} --dna"

    cmd = cmd.format(
        streme=streme_path,
        output=output_dir,
        input_file=kmer_file,
    )
    logging.info("Running STREME: %s", cmd)

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
    with open(raw_gff, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            attributes_dict = {key: value for key, value in (item.split('=') for item in attributes.split(';'))}

            print (f">context_{i}\n{attributes_dict['context']}", file = context_f)
            # print (f">{seqid}_{start}\t{check_existing_motif(attributes_dict)}\n{attributes_dict['context']}", file = context_f)
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
        


### SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L borg
### SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_494_C  Methanoperedens_43_0_circular
#### SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_267_C Methanoperedens_43_25_c

if __name__ == "__main__":
    # prefix = 'break_contigs/break_contigs'
    # prefix = 'seven_contigs/seven_contigs'

    # prefix = 'split_test/split'

    # gff_file = f'{prefix}.reprocess.gff'
    # figure_name = f'{prefix}_heatmap.pdf'
    # matrix_file = f'{prefix}_motif_matrix.csv'
    # context_fasta = f'{prefix}.context.fasta'

    # gff_folder = 'borg/split_bam_dir/'
    # # merge_gff(gff_folder, gff_file)


    prefix = 'test'
    # motif_info_dict, motif_dict = parse_gff(gff_file)
    matrix_file = f'tmp/{prefix}_motif_matrix.csv'
    figure_name = f'tmp/{prefix}_heatmap.pdf'
    motif_info_dict, motif_dict = load_motifs("borg/split_bam_dir/")
    ##motif_info_dict = add_benchmark(motif_info_dict, motif_dict)
    get_motif_matrix(motif_info_dict, motif_dict)
    heatmap(matrix_file)


    # cluster_contigs(matrix_file)
    # run_streme(context_fasta)
    # streme_output = "break_contigs/test.context_streme/"
    # read_streme(streme_output)

    # raw_gff = "/home/shuaiw/methylation/data/borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_10_L.gff"
    # # raw_gff = "/home/shuaiw/methylation/data/borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_11_C.gff"
    # separat_motif(raw_gff)
    # streme_output = raw_gff[:-4] + ".context_streme/"
    # read_streme(streme_output)