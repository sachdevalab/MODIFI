import pandas as pd
import os
from collections import defaultdict
import re
import sys

def read_tsv(anno_file):
    
    df = pd.read_csv(anno_file, sep="\t")
    # Extract contig name from 'Gene' column and add as a new column
    df['contig'] = df['Gene'].apply(lambda x: '_'.join(str(x).split('_')[:-1]))  
    return df

def get_ctg_MTase_motifs(df, bin_name, bins={}):
    ## filter the dataframe for the specific contig
    bin_Mtase = pd.DataFrame()
    if bin_name in bins:
        for ctg in bins[bin_name]:
            df_ctg = df[df['contig'] == ctg]
            # print (df_ctg)
            if not df_ctg.empty:
                bin_Mtase = pd.concat([bin_Mtase, df_ctg])
    else:
        df_ctg = df[df['contig'] == bin_name]
        # print (df_ctg)
        if not df_ctg.empty:
            bin_Mtase = df_ctg
    # print (bin_Mtase)
    return bin_Mtase

def bin_ctgs(bin_file):
    """
    Read the bin file and return a dictionary with contig names as keys and their corresponding bins as values.
    """
    bins = defaultdict(list)
    with open(bin_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                contig_name = parts[0]
                bin_name = parts[1]
                bins[bin_name].append(contig_name)
    return bins

def read_motifs(motif_dir, bin_name, bins={}):
    bin_motif = pd.DataFrame()
    if bin_name in bins:
        for ctg in bins[bin_name]:
            motif_file = os.path.join(motif_dir, f"{ctg}.motifs.csv")
            if os.path.exists(motif_file):
                df_motif = pd.read_csv(motif_file)
                df_motif['contig'] = ctg
                if not df_motif.empty:
                    bin_motif = pd.concat([bin_motif, df_motif])
    else:
        motif_file = os.path.join(motif_dir, f"{bin_name}.motifs.csv")
        if os.path.exists(motif_file):
            df_motif = pd.read_csv(motif_file)
            df_motif['contig'] = bin_name
            if not df_motif.empty:
                bin_motif = df_motif
    return bin_motif

def match_MTase_motifs(bin_Mtase, bin_motif, RM_file):
    ## collect all motifs from bin_motif
    if len(bin_motif) == 0:
        # print(f"No motifs found for {RM_file}. Skipping.")
        return
    motifs = set(bin_motif['motifString'])
    ## enuerate bin_Mtase
    df = pd.DataFrame()
    for index, row in bin_Mtase.iterrows():
        if row['Homolog motif'] in motifs:
            row["match"] = True
            if re.search("Singleton", row["Operon"]):
                row["orphan"] = True
            else:
                row["orphan"] = False
            ## add row to df
            df = pd.concat([df, row.to_frame().T])
    ## output the df to csv
    df.to_csv(RM_file, index=False)

    ## add support info for manual curation
    ## print bin_Mtase and bin_motif to RM_file
    with open(RM_file, 'a') as f:
        f.write("\n# Support info:\n")
        f.write("# bin_Mtase:\n")
        bin_Mtase.to_csv(f, index=False, sep="\t")
        f.write("# bin_motif:\n")
        bin_motif.to_csv(f, index=False, sep="\t")



def RM_main(anno_file, motif_dir, bin_file=None):
    if bin_file:
        bins = bin_ctgs(bin_file)
    else:
        bins = {}
        ## for file in motif_dir:
        for file in os.listdir(motif_dir):
            if file.endswith(".motifs.csv"):
                bin_name = file.split(".motifs.csv")[0]
                bins[bin_name] = [bin_name]

    RM_dir = os.path.join(motif_dir, "..", "RM_system")
    ## create RM_dir if not exists
    if not os.path.exists(RM_dir):
        os.makedirs(RM_dir)

    # bin_name = "B_cereus_971_1"
    for bin_name in bins:
        # print (f"Processing bin: {bin_name}")
        RM_file = os.path.join(RM_dir, f"{bin_name}.RM.csv")
        df = read_tsv(anno_file)
        bin_Mtase = get_ctg_MTase_motifs(df, bin_name, bins)
        bin_motif = read_motifs(motif_dir, bin_name, bins)

        # print (bin_Mtase)
        # print (bin_motif)
        match_MTase_motifs(bin_Mtase, bin_motif, RM_file)

if __name__ == "__main__":
    ## acceept anno_file and motif_dir as arguments
    
    if len(sys.argv) < 3:
        print("Usage: python analyze_RM.py <anno_file> <work_dir>")
        sys.exit(1)
    anno_file = sys.argv[1]
    motif_dir = os.path.join(sys.argv[2], "motifs")
    if len(sys.argv) > 3:
        bin_file = sys.argv[3]
    else:
        bin_file = None

    # anno_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_RM.rm.genes.tsv"
    # motif_dir = "/home/shuaiw/borg/bench/zymo_new_ref2/motifs/"
    RM_main(anno_file, motif_dir, bin_file)

# anno_file = "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_RM.rm.genes.tsv"
# motif_dir = "/home/shuaiw/borg/bench/soil/run1/motifs/"
# bin_name = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_maxbin2_414"
# bin_file = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.bin.tab"


# motif_file = os.path.join(motif_dir, f"{ctg}.motifs.csv")
# bins = bin_ctgs(bin_file)
# df = read_tsv(anno_file)
# bin_Mtase = get_ctg_MTase_motifs(df, bin_name, bins)