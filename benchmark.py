from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt



def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        # return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos
    motif_sites = {}

    for r, contig in REF.items():
        print (r)
        for site in nt_search(str(contig), motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            motif_sites[r + ":" + str(site+exact_pos) + "+"] = motif_new

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            motif_sites[r + ":" + str(site+rev_exact_pos) + "-"] = motif_new

    return motif_sites

def corr_obs_cont1(infer):
    df = pd.read_csv(infer, nrows=1000000)
    df = df[df['coverage'] > 10]
    # df = df[df['kmer_count'] > 10]
    # df = df[df['strand'] == 1]
    ## cal the correlation between tMean and control
    # corr = pearsonr(df['tMean'], df['control'])
    # print ("tMean and control", corr)
    ## plot the distribution of ipd_ratio
    df['ipd_ratio'].plot.hist(bins=100)
    ## svae the plot
    plt.savefig("tmp/ipd_ratio.png")
    ## clean the plt
    # plt.clf()
    ## print the rows with ipd_ratio > 30
    print (df[df['ipd_ratio'] > 10])
    motif_locus_dict = {}
    ## plot the distribution of tMean
    df['tMean'].plot.hist(bins=100)
    ## svae the plot
    plt.savefig("tmp/tMean.png")
    plt.clf()
    ## plot control 
    df['control'].plot.hist(bins=100)
    ## svae the plot
    plt.savefig("tmp/control.png")
    plt.clf()

def corr_obs_cont2(infer):
    df = pd.read_csv(infer)
    # df = df[df['coverage'] > 10]
    modified_ipd_ratio_list = []
    for index, row in df.iterrows():
        if row['strand'] == 1:
            tag = row['refName'] + ":" + str(row['tpl']) + "-"
        else:
            tag = row['refName'] + ":" + str(row['tpl']) + "+"
        if tag in motif_sites:
            modified_ipd_ratio_list.append(row['ipd_ratio'])
    ## plot the distribution of ipd_ratio
    print (len(modified_ipd_ratio_list))
    plt.hist(modified_ipd_ratio_list, bins=100)
    ## svae the plot
    plt.savefig("tmp/modified_ipd_ratio.png")

            

motif_new = "CTGCAG"
exact_pos = 5
my_ref = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
infer = "/home/shuaiw/borg/bench/C227_native/ipd_ratio/CP011331.1.ipd3.csv"
# infer = "/home/shuaiw/borg/bench/ecoli_native/ipd_ratio/CP064387.1.ipd3.csv"

REF = read_ref(my_ref)
# print (REF)
motif_sites = get_motif_sites(REF, motif_new, exact_pos)
# for site in motif_sites:
#     print(site, motif_sites[site])

# print (len(motif_sites))
corr_obs_cont2(infer)