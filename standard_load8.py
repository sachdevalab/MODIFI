"""
load ipd from bam, refer to
https://github.com/PacificBiosciences/kineticsTools/blob/master/kineticsTools/KineticWorker.py
"""


from pbcore.io import AlignmentSet
import numpy as np
import os
import pandas as pd
from collections import defaultdict
# from scipy.stats import norm
import time
import sys
import logging
import pysam

# alignments = None

# Raw ipd record
ipdRec = [('tpl', '<u4'), ('strand', '<i8'), ('ipd', '<f4')]
mapQvThreshold = 0
maxAlignments = 10000 ##1500
randomSeed = None
max_region = 100000


def _subreadNormalizationFactor(rawIpds):
    """
    Normalize subread ipds
    """

    # Default normalization factor -- this value should very rarely get
    # used
    if rawIpds.size < 2:
        return 0.1

    if np.isnan(rawIpds).any():
        print("got nan: %s" % str(rawIpds))

    if rawIpds.mean() < 0.0001:
        print("small")
        print("got small: %s" % str(rawIpds))

    capValue = min(10, np.percentile(rawIpds, 99))
    capIpds = np.minimum(rawIpds, capValue)
    if capIpds.mean() < 0.0001:
        print("small cap")
        print("got small cap: %s" % str(capIpds))
    return capIpds.mean()

def _loadRawIpds(contig_bam, alignments, refGroupId, each_ref, start, end, factor=1.0):
    t0 = time.time()
    # (start, end) = (0, each_ref.Length)

    MIN_IDENTITY = 0.0  
    MIN_READLENGTH = 50
    samfile = pysam.AlignmentFile(contig_bam, "rb", check_sq=False)

    hits = [hit for hit in samfile.fetch("SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C",
                                                        max(start, 0), end)]
    # logging.info("Retrieved %d hits" % len(hits), round(time.time()-t0))
    print ("Retrieved %d hits" % len(hits), "time", round(time.time()-t0))
    if len(hits) > maxAlignments:
        # XXX a bit of a hack - to ensure deterministic behavior when
        # running in parallel, re-seed the RNG before each call
        if randomSeed is None:
            np.random.seed(len(hits))
        hits = np.random.choice(
            hits, size=maxAlignments, replace=False)
    if len(hits) > 0:
        print ("downsample ratio", round(100*maxAlignments/len(hits), 4), "%")
            
    # Maintain separate lists for each strand to speed up sorting
    s0dict = defaultdict(list)
    s1dict = defaultdict(list)
    ipdVect = []

    # for aln in alignments.readsInRange(refGroupId, start, end):
    for aln in hits:
        # Pull out error-free position
        # matched = np.logical_and(np.array(
        #     [x != '-' for x in aln.read()]), np.array([x != '-' for x in aln.reference()]))
        # np.logical_and(np.logical_not(np.isnan(rawIpd)),
        #                 matched, out=matched)
        # normalization = _subreadNormalizationFactor(rawIpd[matched])
        ## check normalization is a valid number or nan
        ## print read name if normalization is nan

        if aln.get_tag("NM") > 3:
            continue
        forward_IPD_info = np.array(aln.get_tag("fi")[::-1])   ## weired, why need to reverse
        reverse_IPD_info = np.array(aln.get_tag("ri"))

        # Get aligned positions on the reference for each read base
        aligned_pairs = aln.get_aligned_pairs(matches_only=True, with_seq=False)
        rawIpd = np.zeros(len(aln.query_sequence))
        matched = np.zeros(len(aln.query_sequence), dtype=bool)
        referencePositions = np.zeros(len(aln.query_sequence), dtype=int)
        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is not None:
                if forward_IPD_info[query_pos] != 0:
                    if ref_pos >= start and ref_pos < end:
                        rawIpd[query_pos] = forward_IPD_info[query_pos]
                        matched[query_pos] = True
                        referencePositions[query_pos] = ref_pos 
            # if not matched[query_pos]:
            #     print ("not matched", aln.query_name, query_pos, ref_pos, start, end, forward_IPD_info[query_pos])
        
        # print ("rawIpd", rawIpd)
        # print ("matched", matched)
        # print ("referencePositions", referencePositions)

        np.logical_and(np.logical_not(np.isnan(rawIpd)),
                        matched, out=matched)
        normalization = _subreadNormalizationFactor(rawIpd[matched])

        if np.isnan(normalization):
            print (f"nan normalization in {aln.query_name}", normalization)
            continue
        if normalization < 0.0001:
            print (f"zero IPD values in {aln.query_name}", normalization, rawIpd[matched])
            continue

        rawIpd /= normalization

        nm = matched.sum()

        # Bail out if we don't have any samples
        if nm == 0:
            continue

        ipd = rawIpd[matched]
        tpl = referencePositions[matched]

        ipdVect += list(ipd)

        for tpl_val, ipd_val in zip(tpl, ipd):
            s0dict[tpl_val].append(ipd_val)
        # else:
        #     for tpl_val, ipd_val in zip(tpl, ipd):
        #         s1dict[tpl_val].append(ipd_val)

    print ("load takes", round(time.time()-t0))
    if len(ipdVect) < 10:
        # Default is there is no coverage
        capValue = 5.0
    else:
        # Compute IPD quantiles on the current block -- will be used for
        # trimming extreme IPDs
        capValue = np.percentile(ipdVect, 99)
    print ("capValue", capValue)
    print ("pos num", len(s0dict), len(s1dict))
    return cal_mean(s0dict, s1dict, each_ref, capValue, start, end, t0)

def cal_mean(s0dict, s1dict, each_ref, capValue, start, end, t0):
    ref_Name = "xx"

    # s0Ipds, s1Ipds = [], []
    result = []
    for pos in range(start, end):
        if pos in s0dict:
            d = np.array(s0dict[pos])
            if len(d) <= 2:
                continue
            coverage = len(d)

            percentile = min(90, (1.0 - 1.0 / (len(d) - 1)) * 100)
            localPercentile = np.percentile(d, percentile)
            local_capValue = max(capValue, localPercentile)

            d = np.minimum(d, local_capValue)

            # Trimmed stats
            tMean = np.mean(d).item()
            tErr = np.std(d).item() / np.sqrt(len(d))
            # print (ref_seq[pos])
            ## 0-based
            result.append([ref_Name, 1, pos, complement_ref_seq[pos], coverage, tMean, tErr])
    
    print ("one strand done", round(time.time()-t0))

    for pos in range(start, end):   
        if pos in s1dict:
            d = np.array(s1dict[pos])
            if len(d) <= 2:
                continue

            coverage = len(d)

            percentile = min(90, (1.0 - 1.0 / (len(d) - 1)) * 100)
            localPercentile = np.percentile(d, percentile)
            local_capValue = max(capValue, localPercentile)
            d = np.minimum(d, local_capValue)

            # Trimmed stats
            tMean = np.mean(d).item()
            tErr = np.std(d).item() / np.sqrt(len(d))
            ## 0-based
            result.append([ref_Name, 0, pos, ref_seq[pos], coverage, tMean, tErr])

    print ("raw ipd is counted", round(time.time()-t0))
    return result

def extract_context(fasta):
    ## load the fasta using biopython
    # print ("loading fasta")
    seq_dict = {}
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        # print(record.id)
        if record.id != "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C":
            continue
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        # seq = str(record.seq)
        ## convert to capital
        seq_dict[record.id] = record.seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
        break
    return seq_dict

def load_IPD(each_ref, contig_bam, df_file):
    # t0 = time.time()
    # alignments = AlignmentSet(contig_bam, referenceFastaFname=fasta)
    # refInfo = alignments.referenceInfoTable
    # # print ("refInfo", refInfo.shape, len(refInfo))
    # for my_ref in refInfo:
    #     if my_ref.Name == each_ref:
    #         each_ref = my_ref
    #         break
    # print ("ref loaded", each_ref.Name, each_ref.Length, round(time.time()-t0))
    # factor = 1.0 / alignments.readGroupTable[0].FrameRate
    # # global alignments
    # refGroupId = alignments.referenceInfo(each_ref.Name).Name

    # # max_length = each_ref.Length
    # max_length = 1000000

    # chunks_num = int(max_length/max_region)
    # chunks_num = 1 if chunks_num < 1 else chunks_num
    
    # print ("chunks_num", chunks_num)
    # result = []
    # for i in range(chunks_num):
    #     start = i*max_region
    #     end = (i+1)*max_region
    #     if end > max_length:
    #         end = max_length
    #     if i == chunks_num-1:
    #         end = max_length
    #     print ("chunk", i, start, end)
    #     result += _loadRawIpds(contig_bam, alignments, refGroupId, each_ref, start, end, factor, )
    #     # break

    # chunk_result = _loadRawIpds(alignments, refGroupId, each_ref, factor, )
    # print ("rawIpds", rawIpds.shape, round(time.time()-t0))
    # combined_df = pd.concat(result, ignore_index=True)
    result = _loadRawIpds(contig_bam, '', '', each_ref, 1, 1000, 1, )
    combined_df = pd.DataFrame(result, columns=['refName', 'strand', 'tpl', 'base', 'coverage', 'tMean', 'tErr'])
    combined_df.to_csv(df_file, index=False)
    print ("raw ipd df saved", df_file, round(time.time()))






if __name__ == "__main__":

    # subread_bam = "/home/shuaiw/methylation/data/borg/b_contigs/11.align.bam"
    # fasta = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa"
    # subread_bam = "/home/shuaiw/methylation/data/borg/new_test4/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_267_C.bam"
    # fasta = "/home/shuaiw/methylation/data/borg/new_test4/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_267_C.fa"
    # outputfile = "/home/shuaiw/methylation/data/borg/b_contigs/test/test.csv"

    # subread_bam = "/home/shuaiw/methylation/data/borg/all_borg/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.filter.bam"
    # fasta = "/home/shuaiw/methylation/data/borg/all_borg/all_borg.fasta"

    # subread_bam = "/home/shuaiw/methylation/data/borg/human/human_191315.ccs.align.bam"
    # fasta = "/home/shuaiw/borg/hg38/GCF_000001405.40_GRCh38.p14_genomic.fasta"

    subread_bam = "/home/shuaiw/methylation/data/borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam"
    fasta = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
    outputfile =  "/home/shuaiw/methylation/data/borg/all_borg/test_ccs.csv"

    # subread_bam = sys.argv[1]
    # fasta = sys.argv[2]
    # outputfile = sys.argv[3]

    # ## set default value for maxAlignments if not set
    # if len(sys.argv) > 4:
    #     maxAlignments = int(sys.argv[4])
    #     print ("para maxAlignments", maxAlignments)

    ## build outdir if not exists
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)

    seq_dict = extract_context(fasta)
    print ("fasta loaded, contig num:", len(seq_dict))

    for each_ref in seq_dict:
        ref_seq = seq_dict[each_ref]
        ## complement the sequence
        complement_ref_seq = ref_seq.complement()
        load_IPD(each_ref, subread_bam, outputfile)

    # print ("IPD loaded")

    # python /home/shuaiw/Methy/standard_load6.py /home/shuaiw/methylation/data/borg/new_test5/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.bam /home/shuaiw/methylation/data/borg/new_test5/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.fa /home/shuaiw/methylation/data/borg/new_test5/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.ipd1.csv

    ## plot
    data = []
    df = pd.read_csv(outputfile)
    for index, row in df.iterrows():
        if row['tpl'] < 200:
            data.append([row['tpl'], row['base'], row['tMean'], 1-row['strand'], 'hifi'])
    
    df2 = pd.read_csv("/home/shuaiw/methylation/data/borg/all_test/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.ipd1.csv")
    df2 = df2[df2['tpl'] < 200]
    for index, row in df2.iterrows():
        data.append([row['tpl'], row['base'], row['tMean'], row['strand'], 'subread'])
    
    ### plot line plot
    df = pd.DataFrame(data, columns = ["tpl", "base", "tMean", "strand", "type"])
    import seaborn as sns
    import matplotlib.pyplot as plt
    sns.set(style="whitegrid")
    fig, axs = plt.subplots(2, 1, figsize=(20, 10))  # Create a single row with two columns
    ## plot for each strand separately
    sns.lineplot(data=df[df['strand'] == 1], x="tpl", y="tMean", hue="type", style="strand", ax=axs[0])
    sns.lineplot(data=df[df['strand'] == 0], x="tpl", y="tMean", hue="type", style="strand", ax=axs[1])
    ## save the plot in tmp
    plt.savefig("tmp/ipd_lineplot.png")
