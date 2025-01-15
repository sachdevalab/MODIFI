import pysam
import numpy as np
from collections import defaultdict
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import random
import pickle
import os
# import pbcore.io as pb
# from scipy.stats import ttest_ind

def cappedSlog(v, exclude=99):
    q = np.percentile(v, exclude)
    v2 = v.copy()
    # v2 = v2[~np.isnan(v2)]
    v2[v2 > q] = q
    v2[v2 <= 0] = 1. / (75 + 1)
    # return np.log(v2)
    return v2

def cappedSlog2(v, exclude=99):
    ## max value is 95% percentile
    q = np.percentile(v, exclude)
    v2 = v.copy()
    ## replace the value larger than 95% percentile with 95% percentile
    print (q)
    v2[v2 > q] = q
    return v2

def compare_IPD(forward_IPD_info, reverse_IPD_info):
    """
    diff 141207 125951
    diff 239076 219284
    diff 199569 180695
    diff 163348 147796
    diff 222323 191897
    diff 161598 139454
    diff 144831 131893
    diff 253007 230459
    diff 170352 162346
    forward_IPD_info and reverse_IPD_info are the IPD information for the same read, but in different directions.
    """
    diff1 = 0
    diff2 = 0
    for i in range(len(forward_IPD_info)):
        diff1 += abs(forward_IPD_info[i] - reverse_IPD_info[i])
        diff2 += abs(forward_IPD_info[i] - reverse_IPD_info[-i])
    print ("diff", diff1, diff2)

def parse_gff(file_path):
    benchmark_modified_pos = {}
    motif_info_dict = defaultdict(dict)
    test_modified_pos = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            attributes_dict = {key: value for key, value in (item.split('=') for item in attributes.split(';'))}
            score = int(score)
            if score > 100:
                test_modified_pos.append([seqid, int(start), score, attributes_dict['IPDRatio'], feature_type, strand])
            # if feature_type != 'm6A':
            #     continue
            if score > 0:
                benchmark_modified_pos[int(start)] = [seqid, int(start), score, attributes_dict['IPDRatio'], feature_type, strand]
            # if len(test_modified_pos) > 10:
            #     break
    return test_modified_pos, benchmark_modified_pos

def read_subread_bam(bam_file):
    """
    see this for fi and ri tag : 
    https://pacbiofileformats.readthedocs.io/en/12.0/BAM.html#use-of-read-tags-for-hifi-per-read-base-kinetic-information
    """
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    f = open(f'/home/shuaiw/methylation/data/borg/human//human_observed_ipd.csv', 'w')
    print ("start reading")

    ### for each contig in the bam file

    for contig, length in zip(samfile.references, samfile.lengths):
        print(f"Processing contig: {contig}, Length: {length}")


        contig_forward_dict = defaultdict(list)
        contig_reverse_dict = defaultdict(list)

        i = 0
        for read in samfile.fetch(contig):

            IPD_info = read.get_tag("ip")

            # Get aligned positions on the reference for each read base
            aligned_pairs = read.get_aligned_pairs(matches_only=True, with_seq=False)
            for query_pos, ref_pos in aligned_pairs:
                if ref_pos is not None:
                    if IPD_info[query_pos] != 0:
                        
                        if read.is_reverse == False:
                            contig_reverse_dict[ref_pos].append(IPD_info[query_pos])
                        else:
                            ## get the reverse index for query position
                            query_pos = len(IPD_info) - query_pos - 1
                            contig_forward_dict[ref_pos].append(IPD_info[query_pos])
            i += 1
            if i % 100000000 == 0:
                print (i)

        observed_IPD_list = get_IPD_list(length, contig_forward_dict)
        observed_IPD_reverse_list = get_IPD_list(length, contig_reverse_dict)
        
        for i in range(length):
            f.write(f'{contig},{i},{observed_IPD_list[i]},{observed_IPD_reverse_list[i]}\n')
        break
    f.close()

def read_hifi_bam(bam_file):
    """
    see this for fi and ri tag : 
    https://pacbiofileformats.readthedocs.io/en/12.0/BAM.html#use-of-read-tags-for-hifi-per-read-base-kinetic-information
    """
    samfile = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    contig_forward_dict = defaultdict(list)
    contig_reverse_dict = defaultdict(list)
    for read in samfile.fetch("SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C", 1, 1000):
        ## if NM tag larger than 3, skip
        if read.get_tag("NM") > 3:
            continue

        forward_IPD_info = read.get_tag("fi")[::-1]   ## weired, why need to reverse
        reverse_IPD_info = read.get_tag("ri")

        # Get aligned positions on the reference for each read base
        aligned_pairs = read.get_aligned_pairs(matches_only=True, with_seq=False)
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is not None:
                if forward_IPD_info[query_pos] != 0:
                    contig_forward_dict[ref_pos].append(forward_IPD_info[query_pos])
                if reverse_IPD_info[query_pos] != 0:
                    contig_reverse_dict[ref_pos].append(reverse_IPD_info[query_pos])

    observed_IPD_list = get_IPD_list(1000, contig_forward_dict)
    observed_IPD_reverse_list = get_IPD_list(1000, contig_reverse_dict)
    return observed_IPD_list, observed_IPD_reverse_list

def read_bam_pb(bam_file):
    samfile = pb.IndexedBamReader(bam_file)
    ## enumerate the reads on the ref
    for read in samfile.fetch("SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C", 1, 1000):
        ## if NM tag larger than 3, skip
        if read.get_tag("NM") > 3:
            continue

def baseline_correction(IPD_list, reverse_IPD_list):
    # IPD_list = cappedSlog2(np.array(IPD_list))
    # reverse_IPD_list = cappedSlog2(np.array(reverse_IPD_list))

    # overall_mean = np.median(list(IPD_list)+list(reverse_IPD_list))
    overall_mean = np.mean(list(IPD_list)+list(reverse_IPD_list))
    IPD_list = np.array(IPD_list)
    reverse_IPD_list = np.array(reverse_IPD_list)



    IPD_list = IPD_list/overall_mean
    reverse_IPD_list = reverse_IPD_list/overall_mean

    # IPD_list = IPD_list/np.mean(IPD_list)
    # reverse_IPD_list = reverse_IPD_list/np.mean(reverse_IPD_list)

    return IPD_list, reverse_IPD_list


def get_IPD_list(contig_length, contig_dict):
    average_IPD = []
    for pos in range(contig_length):
        if pos not in contig_dict:
            average_IPD.append(0)
        else:
            average_IPD.append(round(np.mean(contig_dict[pos]), 1))
    return average_IPD

class Contig():

    def __init__(self):
        self.contig_name = None
        self.contig_length = None
        self.contig_dict = None
        self.average_IPD = None
        self.average_IPD_reverse = None
        self.all_IPD = None
        self.all_IPD_reverse = None
    
    def get_average_IPD(self, contig_name, contig_length, contig_dict, contig_reverse_dict):
        self.contig_name = contig_name
        self.contig_length = contig_length
        self.contig_dict = contig_dict
        self.average_IPD = []
        self.average_IPD_reverse = []
        self.all_IPD= []
        self.all_IPD_reverse = []
        for pos in range(self.contig_length):
            if pos not in self.contig_dict:
                self.average_IPD.append(0)
                self.average_IPD_reverse.append(0)
                self.all_IPD.append([0])
                self.all_IPD_reverse.append([0])
            else:
                self.average_IPD.append(round(np.mean(self.contig_dict[pos]), 1))
                self.all_IPD.append(self.contig_dict[pos])
                self.average_IPD_reverse.append(round(np.mean(contig_reverse_dict[pos]), 1))
                self.all_IPD_reverse.append(contig_reverse_dict[pos])
                if pos > 430 and pos < 460:
                    print(pos, round(np.mean(self.contig_dict[pos])), round(np.mean(contig_reverse_dict[pos]), 1))

    def play_ipd(self):
        # print(self.average_IPD[:50])
        ## seqid, int(start), score, attributes_dict['IPDRatio']
        data = []
        ## shuffle the test_modified_pos
        random.shuffle(test_modified_pos)
        j = 0
        for modified_pos in test_modified_pos:
            
            for i in range(modified_pos[1]-40, modified_pos[1]+40):
                data.append([i, self.average_IPD[i], modified_pos[1], modified_pos[4], '+'])
                data.append([i, self.average_IPD_reverse[i], modified_pos[1], modified_pos[4], '-'])
            j += 1
            if j > 4:
                break
        ## transform data to dataframe
        df = pd.DataFrame(data, columns=['pos', 'IPD', 'modified', 'feature_type', 'strand'])
        ## output the dataframe to csv
        df.to_csv(f'customized/df.csv', index=False)
        ## plot line plot with seaborn
        # Plot the line plot with seaborn
        # plt.figure(figsize=(12, 6))
        # sns.lineplot(x='pos', y='IPD', data=df)
        # plt.title('IPD Values Along the Reference Genome')
        # plt.xlabel('Reference Position')
        # plt.ylabel('IPD')
        # ## save the plot
        # plt.savefig(f'line_{modified_pos[1]}.pdf')
        # break

    def play_ipd2(self):
        # print(self.average_IPD[:50])
        ## seqid, int(start), score, attributes_dict['IPDRatio']
        data = []

        start = 400
        end = 600
        for i in range(start, end):
            if i + 1 in benchmark_modified_pos:
                truth_modified = benchmark_modified_pos[i + 1]
                strand = truth_modified[5]
                truth_type = truth_modified[4]
                score = truth_modified[2]
                if strand == '+':
                    data.append([i, self.average_IPD[i], truth_type, '+', score])
                    data.append([i, self.average_IPD_reverse[i], 'NA', '-', 1])
                else:
                    data.append([i, self.average_IPD_reverse[i], truth_type, '-', score])
                    data.append([i, self.average_IPD[i], 'NA', '+', 1])
            else:
                data.append([i, self.average_IPD[i], 'NA', '+', 1])
                data.append([i, self.average_IPD_reverse[i], 'NA', '-', 1])

        ## transform data to dataframe
        df = pd.DataFrame(data, columns=['pos', 'IPD', 'feature_type', 'strand', 'score'])
        ## output the dataframe to csv
        df.to_csv(f'customized/bin_df.csv', index=False)
        ## plot line plot with seaborn

    def play_ipd3(self, IPD_ratio_list, reverse_IPD_ratio_list, control_list, \
                  reverse_control_list, observed_IPD_list, observed_IPD_reverse_list,\
                summary_IPD_list, summary_reverse_IPD_list):

        data = []

        start = 400
        end = 700
        for i in range(start, end):
            forward_type, reverse_type = 'NA', 'NA'
            forward_score, reverse_score = 1, 1
            if i + 1 in benchmark_modified_pos:
                truth_modified = benchmark_modified_pos[i + 1]
                strand = truth_modified[5]
                truth_type = truth_modified[4]
                score = truth_modified[2]
                if strand == '+':
                    forward_type = truth_type
                    forward_score = score
                else:
                    reverse_type = truth_type
                    reverse_score = score
            data.append([i, summary_IPD_list[i], forward_type, '+ IPD_summary IPD', forward_score])
            data.append([i, summary_reverse_IPD_list[i], reverse_type, '- IPD_summary IPD', reverse_score])
            data.append([i, IPD_ratio_list[i], forward_type, '+ IPD_Ratio', forward_score])
            data.append([i, reverse_IPD_ratio_list[i], reverse_type, '- IPD_Ratio', reverse_score])
            data.append([i, control_list[i], forward_type, '+ Control', forward_score])
            data.append([i, reverse_control_list[i], reverse_type, '- Control', reverse_score])
            data.append([i, observed_IPD_list[i], forward_type, '+ Our Observed', forward_score])
            data.append([i, observed_IPD_reverse_list[i], reverse_type, '- Our Observed', reverse_score])


        ## transform data to dataframe
        df = pd.DataFrame(data, columns=['pos', 'IPD', 'feature_type', 'strand', 'score'])
        ## output the dataframe to csv
        df.to_csv(line_data, index=False)
        ## plot line plot with seaborn
    
    def plot_ipd_dist(self):
        ## plot the distribution of IPD values
        plt.figure(figsize=(12, 6))
        sns.histplot(self.average_IPD, bins=100)
        plt.title('Distribution of IPD Values')
        plt.xlabel('IPD')
        plt.ylabel('Frequency')
        plt.savefig('IPD_dist.pdf')

    def simple_cutoff(self):
        ## simple cutoff method
        inferred_modified = {}
        for pos in range(len(self.average_IPD)):
            if self.average_IPD[pos] > 40:
                inferred_modified[pos + 1] = self.average_IPD[pos]   ## convert 0-based to 1-based

        for pos in range(len(self.average_IPD_reverse)):
            if self.average_IPD_reverse[pos] > 40:
                inferred_modified[pos + 1] = self.average_IPD_reverse[pos]   ## convert 0-based to 1-based
        ## examine the inferred modified positions
        true_identified = 0
        false_identified = 0
        for infer in inferred_modified:
            if infer in benchmark_modified_pos:
                true_identified += 1
            else:
                false_identified += 1
        missed = 0
        for true in benchmark_modified_pos:
            if true not in inferred_modified:
                print ('missed', true, self.average_IPD[true-1])
                missed += 1
        print (len(inferred_modified), len(benchmark_modified_pos), "true_identified", true_identified,\
                "false_identified", false_identified, "missed", missed)
        print ("true positive rate", true_identified/len(benchmark_modified_pos))
        print ("false positive rate", false_identified/len(inferred_modified))
        print ("missed rate", missed/len(benchmark_modified_pos))
    
    def store_IPD(self):
        ## use pickle to store the IPD information
        
        with open('customized/IPD.pkl', 'wb') as f:
            pickle.dump(self.average_IPD, f)

    def load_IPD(self):
        ## load the IPD information
        with open('customized/IPD.pkl', 'rb') as f:
            self.average_IPD = pickle.load(f)
                
def extract_context(fasta):
    ## load the fasta using biopython
    print ("loading fasta")
    seq_dict = {}
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        # print(record.id)
        if record.id != "NC_000001.11":
            continue
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)
        ## convert to capital
        raw_seq = seq.upper()
        seq_dict[record.id] = raw_seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
    return seq_dict

def prepare_data(seq, control_list, up=8, down=4):
    # print(seq[:100])
    # seq, control_list = seq[:1000], control_list[:1000]
    kmer_baseline_dict = defaultdict(list)
    # y = control_list[up:len(seq) - down]
    for i in range(up, len(seq) - down):
        kmer = seq[i-up:i+down]
        ipd = control_list[i]
        kmer_baseline_dict[kmer].append(ipd)
    return kmer_baseline_dict



def extract_ipd_ratio(file_path):
    ## open it using pandas
    IPD_ratio_list, reverse_IPD_ratio_list = [], []
    IPD_list, reverse_IPD_list = [], []
    control_list, reverse_control_list = [], []
    df = pd.read_csv(file_path)
    ## extract the IPD ratio
    # i = 0
    # for index, row in df.iterrows():
    #     if int(row['strand']) == 0:
    #         IPD_ratio_list.append(float(row['ipdRatio']))
    #         IPD_list.append(float(row['tMean']))
    #         control_list.append(float(row['modelPrediction']))
    #     elif int(row['strand']) == 1:
    #         reverse_IPD_ratio_list.append(float(row['ipdRatio']))
    #         reverse_IPD_list.append(float(row['tMean']))
    #         reverse_control_list.append(float(row['modelPrediction']))
    #     else:
    #         print (row['strand'], "error", int(row['strand']), int(row['strand']) == 0)
    #     i += 1
    #     if i > 2000:
    #         break
    # return IPD_ratio_list, reverse_IPD_ratio_list, IPD_list, reverse_IPD_list, control_list, reverse_control_list
    return df.loc[df['strand'] == 0, 'ipdRatio'].tolist(), df.loc[df['strand'] == 1, 'ipdRatio'].tolist(),\
            df.loc[df['strand'] == 0, 'tMean'].tolist(), df.loc[df['strand'] == 1, 'tMean'].tolist(),\
            df.loc[df['strand'] == 0, 'modelPrediction'].tolist(), df.loc[df['strand'] == 1, 'modelPrediction'].tolist()


subread_bam = "/home/shuaiw/methylation/data/borg/human/human_000733.subreads.align.bam"
read_subread_bam(subread_bam)






