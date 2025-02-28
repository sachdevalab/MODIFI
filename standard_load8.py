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


def cal_mean(s0dict, s1dict, each_ref, capValue, start, end, t0):
    ref_Name = each_ref

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

def _loadRawIpds_hifi(contig_bam, alignments, refGroupId, each_ref, start, end, factor=1.0):
    t0 = time.time()
    # (start, end) = (0, each_ref.Length)

    MIN_IDENTITY = 0.0  
    MIN_READLENGTH = 50
    # samfile = pysam.AlignmentFile(contig_bam, "rb", check_sq=False)
    # hits = [hit for hit in samfile.fetch("SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C",
    #                                                     max(start, 0), end)]

    hits = [hit for hit in alignments.fetch(each_ref,
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
        if aln.get_tag("NM") > 3:
            continue
        forward_IPD_info = np.array(aln.get_tag("fi")[::-1]) * factor  ## weired, why need to reverse
        reverse_IPD_info = np.array(aln.get_tag("ri")) * factor

        # Initialize arrays
        rawIpd = np.zeros(len(aln.query_sequence))
        matched = np.zeros(len(aln.query_sequence), dtype=bool)
        referencePositions = np.zeros(len(aln.query_sequence), dtype=int)

        rev_rawIpd = np.zeros(len(aln.query_sequence))
        rev_matched = np.zeros(len(aln.query_sequence), dtype=bool)
        rev_referencePositions = np.zeros(len(aln.query_sequence), dtype=int)

        """
        # Get aligned positions on the reference for each read base
        aligned_pairs = aln.get_aligned_pairs(matches_only=True, with_seq=False)

        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is not None:
                if forward_IPD_info[query_pos] != 0:
                    if ref_pos >= start and ref_pos < end:
                        rawIpd[query_pos] = forward_IPD_info[query_pos]
                        matched[query_pos] = True
                        referencePositions[query_pos] = ref_pos 
                if reverse_IPD_info[query_pos] != 0:
                    if ref_pos >= start and ref_pos < end:
                        rev_rawIpd[query_pos] = reverse_IPD_info[query_pos]
                        rev_matched[query_pos] = True
                        rev_referencePositions[query_pos] = ref_pos
        """

        # Get aligned positions on the reference for each read base
        aligned_pairs = np.array(aln.get_aligned_pairs(matches_only=True, with_seq=False))

        # Extract query and reference positions as NumPy arrays
        query_positions = aligned_pairs[:, 0]
        reference_positions = aligned_pairs[:, 1]

        # Mask for valid reference positions within the start-end range
        valid_mask = (reference_positions >= start) & (reference_positions < end)

        # Apply vectorized filtering
        forward_mask = valid_mask & (forward_IPD_info[query_positions] != 0)
        reverse_mask = valid_mask & (reverse_IPD_info[query_positions] != 0)

        # Assign values efficiently
        rawIpd[query_positions[forward_mask]] = forward_IPD_info[query_positions[forward_mask]]
        matched[query_positions[forward_mask]] = True
        referencePositions[query_positions[forward_mask]] = reference_positions[forward_mask]

        rev_rawIpd[query_positions[reverse_mask]] = reverse_IPD_info[query_positions[reverse_mask]]
        rev_matched[query_positions[reverse_mask]] = True
        rev_referencePositions[query_positions[reverse_mask]] = reference_positions[reverse_mask]



        ipd, tpl = norm(rawIpd, referencePositions, matched, aln)
        rev_ipd, rev_tpl = norm(rev_rawIpd, rev_referencePositions, rev_matched, aln)
        if len(ipd) == 0 or len(rev_ipd) == 0:
            continue

        ipdVect += list(ipd)
        ipdVect += list(rev_ipd)

        for tpl_val, ipd_val in zip(tpl, ipd):
            s1dict[tpl_val].append(ipd_val)

        for tpl_val, ipd_val in zip(rev_tpl, rev_ipd):
            s0dict[tpl_val].append(ipd_val)

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

def norm(rawIpd, referencePositions, matched, aln):
    np.logical_and(np.logical_not(np.isnan(rawIpd)),
                    matched, out=matched)
    normalization = _subreadNormalizationFactor(rawIpd[matched])

    if np.isnan(normalization):
        print (f"nan normalization in {aln.query_name}", normalization)
        return [], []
    if normalization < 0.0001:
        print (f"zero IPD values in {aln.query_name}", normalization, rawIpd[matched])
        return [], []

    rawIpd /= normalization

    nm = matched.sum()

    # Bail out if we don't have any samples
    if nm == 0:
        return [], []

    ipd = rawIpd[matched]
    tpl = referencePositions[matched]
    return ipd, tpl


def load_IPD_hifi(each_ref, ref_seq, contig_bam, df_file):
    print (f"handle contig {each_ref}, with length {len(ref_seq)}...")
    t0 = time.time()
    alignments = pysam.AlignmentFile(contig_bam, "rb", check_sq=False)
    read_groups = alignments.header.get("RG", [])
    if read_groups:
        # Get Frame Rate from the first Read Group (if present)
        frame_rate = float(read_groups[0].get("FR", 1.0))  # Default to 1.0 if not found
        factor = 1.0 / frame_rate  # Calculate normalization factor
        print(f"Frame Rate: {frame_rate}")
        print(f"Factor: {factor}")
    else:
        factor = 1.0

    # # max_length = each_ref.Length
    # max_length = 100000
    max_length = len(ref_seq)
    refGroupId = each_ref

    chunks_num = int(max_length/max_region)
    chunks_num = 1 if chunks_num < 1 else chunks_num
    
    print ("chunks_num", chunks_num)
    result = []
    for i in range(chunks_num):
        start = i*max_region
        end = (i+1)*max_region
        if end > max_length:
            end = max_length
        if i == chunks_num-1:
            end = max_length
        print ("chunk", i, start, end)
        result += _loadRawIpds_hifi(contig_bam, alignments, refGroupId, each_ref, start, end, factor, )
        # break
    combined_df = pd.DataFrame(result, columns=['refName', 'strand', 'tpl', 'base', 'coverage', 'tMean', 'tErr'])
    combined_df.to_csv(df_file, index=False)
    print ("raw ipd df saved", df_file, round(time.time()-t0))






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
        load_IPD_hifi(each_ref, seq_dict[each_ref], subread_bam, outputfile)

    # print ("IPD loaded")

    # python /home/shuaiw/Methy/standard_load6.py /home/shuaiw/methylation/data/borg/new_test5/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.bam /home/shuaiw/methylation/data/borg/new_test5/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.fa /home/shuaiw/methylation/data/borg/new_test5/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.ipd1.csv

    ## plot
    data = []
    df = pd.read_csv(outputfile)
    for index, row in df.iterrows():
        if row['tpl'] < 300 and row['tpl'] > 200:
            data.append([row['tpl'], row['base'], row['tMean'], row['strand'], 'hifi'])
    
    df2 = pd.read_csv("/home/shuaiw/methylation/data/borg/all_test/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_105_C.ipd1.csv")
    df2 = df2[(df2['tpl'] < 300) & (df2['tpl'] > 200)]
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
