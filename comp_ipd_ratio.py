import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict
import numpy as np
# from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import sys
import seaborn as sns
import matplotlib.pyplot as plt



def p_value_right_tail(x, mu, sigma):
    """
    Calculate the one-tailed p-value for x being in the right tail of a normal distribution.

    Parameters:
    x (float): The observed value
    mu (float): The mean of the normal distribution
    sigma (float): The standard deviation of the normal distribution

    Returns:
    float: The p-value for x being in the right tail
    """
    # Calculate the z-score
    z = (x - mu) / sigma
    
    # Calculate the right-tail p-value
    p_value = 1 - norm.cdf(z)
    
    return p_value

def calculate_x_from_pvalue(p_value, mu, sigma, tail="right"):
    """
    Calculate the value(s) x that correspond to a given p-value in a normal distribution.

    Parameters:
    p_value (float): The p-value (probability)
    mu (float): The mean of the normal distribution
    sigma (float): The standard deviation of the normal distribution
    tail (str): The tail type. Options are "right", "left", or "two-sided".

    Returns:
    float or tuple: The value(s) x corresponding to the given p-value
    """
    if tail == "right":
        # Right-tail: Find the value where P(Z >= z) = p_value
        z = norm.ppf(1 - p_value)
        x = mu + z * sigma
        return x

def get_ipd_ratio(csv, output, gff, figure_file, min_cov=5):
    df = pd.read_csv(csv, sep = ",")
    df['tpl'] = df['tpl'].astype(int)
    df['strand'] = df['strand'].astype(int)
    df['coverage'] = df['coverage'].astype(int)
    df['kmer_count'] = df['kmer_count'].astype(int)

    ## remove the elements with ipd_ratio == 0
    df = df[df['ipd_ratio'] != 0]
    ## remove the elements with coverage < min_cov
    df = df[df['coverage'] >= int(min_cov)]

    mean = df['ipd_ratio'].mean()
    std = df['ipd_ratio'].std()
    df['pvalue'] = df['ipd_ratio'].apply(lambda x: p_value_right_tail(x, mean, std))

    ## add score column, score i s -10logpvalue
    ## a Phred-transformed QV, QV =−10 log10 p
    # df['score'] = -10 * np.log10(df['pvalue'])

    ### round all the float values to 2 decimal places

    df = df.round(4)

    df.to_csv(output, index=False)
    visu(df, figure_file)
    ## ImportError: Matplotlib requires numpy>=1.23; you have 1.22.4
    get_gff(df, gff)

def visu(df, figure_path):
    if df.empty:
        print("DataFrame is empty. Cannot generate visualization.")
        ## construct an  empty figure
        plt.figure()
        plt.savefig(figure_path)
        return
    
    sns.set(style="whitegrid")
    fig, axs = plt.subplots(2, 2, figsize=(20, 10))
    ## first row is covergae, second is tMean, third is control, fourth is ipd_ratio
    ## plot for each strand separately
    sns.histplot(df, x="coverage", hue="strand", multiple="stack", ax=axs[0, 0])
    sns.histplot(df, x="tMean", hue="strand", multiple="stack", ax=axs[0, 1])
    sns.histplot(df, x="control", hue="strand", multiple="stack", ax=axs[1, 0])
    sns.histplot(df, x="ipd_ratio", hue="strand", multiple="stack", ax=axs[1, 1])
    ## save the plot

    ### cal the number with strand 1 and 0
    # print (df[df['strand'] == 1].shape[0])
    # print (df[df['strand'] == 0].shape[0])
    plt.savefig(figure_path)

def phred_qv(p, max_qv=60):
    """Compute Phred-transformed Quality Value, handling p=0 cases."""
    if p == 0:
        return max_qv  # Cap at a defined maximum QV
    return min(round(-10 * np.log10(p)), max_qv)

def get_ref(ref):
    # print ("loading fasta")
    seq_dict = {}
    from Bio import SeqIO
    for record in SeqIO.parse(ref, "fasta"):
        seq_dict[record.id] = record.seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
    return seq_dict


def get_gff(df, gff, p_cut = 0.05):
    print ("start get gff...")
    df = df[df['pvalue'] <= p_cut]
    # gff = "tmp/test.gff"
    f = open(gff, "w")
    print ("##gff-version 3", file = f)
    print ("##source-version our_method", file = f)
    print ("##source-commandline xxx", file = f)
    for seq in seq_dict:
        print ("##sequence-region\t" + seq + "\t1\t" + str(len(seq_dict[seq])), file = f)
    # print ("##sequence-region\tCP064388.1\t1\t75554", file = f )
    for index, row in df.iterrows():
        score = phred_qv(row['pvalue'])
        one_based_tpl = row['tpl'] + 1
        if row['strand'] == 1:
            my_strand = '-'
        else:
            my_strand = '+'
        # Handle context sequence with padding
        start = max(0, row['tpl'] - 20)
        end = min(len(seq_dict[row['refName']]), row['tpl'] + 21)
        context = seq_dict[row['refName']][start:end]
        if start == 0:
            context = 'N' * (20 - row['tpl']) + context
        if end == len(seq_dict[row['refName']]):
            context = context + 'N' * (row['tpl'] + 21 - len(seq_dict[row['refName']]))
        
        if row['strand'] == 1:
            context = context.reverse_complement()
        
        anno = f"coverage={row['coverage']};context={context};IPDRatio={row['ipd_ratio']}"
        print (f"{row['refName']}\tkinModCall\tmodified_base\t{one_based_tpl}\t{one_based_tpl}\t{score}\t{my_strand}\t.\t{anno}", file = f)
    f.close()


if __name__ == "__main__":
    csv = sys.argv[1]
    output = sys.argv[2]
    gff = sys.argv[3]
    # ref = "/home/shuaiw/borg/bench/ecoli_native/contigs/CP064388.1.fa"
    ref = sys.argv[4]
    # fig = None #sys.argv[5]
    figure_file = sys.argv[5]
    min_cov = sys.argv[6]
    seq_dict = get_ref(ref)
    print ("loaded fasta")
    get_ipd_ratio(csv, output, gff, figure_file, min_cov)


    # csv = "/home/shuaiw/methylation/data/borg/b_contigs/control/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_0_219069.ipd2.csv"
    # output = "/home/shuaiw/methylation/data/borg/b_contigs/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_0_219069.ipd3.csv"
    # get_ipd_ratio(csv, output)
    # get_raw_ipd(df)

    # python comp_ipd_ratio.py /home/shuaiw/borg/bench/ecoli_native/control/CP064388.1.ipd2.csv tmp/test.csv tmp/test.gff /home/shuaiw/borg/bench/ecoli_native/contigs/CP064388.1.fa
    # ~/smrtlink/motifMaker find -f /home/shuaiw/borg/bench/ecoli_native/contigs/CP064388.1.fa -g tmp/test.gff -j 1 -o tmp/test.motif.csv -m 0
