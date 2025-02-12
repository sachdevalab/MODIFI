"""
load ipd from bam, refer to
https://github.com/PacificBiosciences/kineticsTools/blob/master/kineticsTools/KineticWorker.py
"""


from pbcore.io import AlignmentSet
import numpy as np
import os
import pandas as pd
from collections import defaultdict
import pickle
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import time
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
from multiprocessing import Manager, Process, Lock
import sys
import logging

# alignments = None

# Raw ipd record
ipdRec = [('tpl', '<u4'), ('strand', '<i8'), ('ipd', '<f4')]
mapQvThreshold = 0
maxAlignments = 1500
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
    return capIpds.mean()

def _loadRawIpds(alignments, refGroupId, each_ref, start, end, factor=1.0):
    t0 = time.time()
    # (start, end) = (0, each_ref.Length)

    MIN_IDENTITY = 0.0  
    MIN_READLENGTH = 50

    hits = [hit for hit in alignments.readsInRange(refGroupId,
                                                        max(start, 0), end)
            if ((hit.mapQV >= mapQvThreshold) and
                (hit.identity >= MIN_IDENTITY) and
                (hit.readLength >= MIN_READLENGTH))]
    # logging.info("Retrieved %d hits" % len(hits), time.time()-t0)
    print ("Retrieved %d hits" % len(hits), time.time()-t0)
    if len(hits) > maxAlignments:
        # XXX a bit of a hack - to ensure deterministic behavior when
        # running in parallel, re-seed the RNG before each call
        if randomSeed is None:
            np.random.seed(len(hits))
        hits = np.random.choice(
            hits, size=maxAlignments, replace=False)
            
    # Maintain separate lists for each strand to speed up sorting
    s0dict = defaultdict(list)
    s1dict = defaultdict(list)
    ipdVect = []

    # for aln in alignments.readsInRange(refGroupId, start, end):
    for aln in hits:
        # Pull out error-free position
        matched = np.logical_and(np.array(
            [x != '-' for x in aln.read()]), np.array([x != '-' for x in aln.reference()]))

        # Normalize kinetics of the entire subread
        rawIpd = aln.IPD() * factor

        np.logical_and(np.logical_not(np.isnan(rawIpd)),
                        matched, out=matched)

        normalization = _subreadNormalizationFactor(rawIpd[matched])
        rawIpd /= normalization

        # Trim down to just the position that cover our interval
        referencePositions = aln.referencePositions()
        np.logical_and(referencePositions < end, matched, matched)
        np.logical_and(referencePositions >= start, matched, matched)
        nm = matched.sum()

        # Bail out if we don't have any samples
        if nm == 0:
            continue

        ipd = rawIpd[matched]
        tpl = referencePositions[matched]

        ipdVect += list(ipd)

        if aln.isForwardStrand:
            for tpl_val, ipd_val in zip(tpl, ipd):
                s0dict[tpl_val].append(ipd_val)
        else:
            for tpl_val, ipd_val in zip(tpl, ipd):
                s1dict[tpl_val].append(ipd_val)

    print ("load takes", time.time()-t0)
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
    ref_Name = each_ref.Name

    # s0Ipds, s1Ipds = [], []
    result = []
    for pos in range(start, end):
        if pos in s0dict:
            d = np.array(s0dict[pos])
            if len(d) <= 2:
                continue

            # NOTE -- this is where the strand flipping occurs -- make sure to
            # reproduce this in the all calling methods
            strand = 1 - 0
            coverage = len(d)

            percentile = min(90, (1.0 - 1.0 / (len(d) - 1)) * 100)
            localPercentile = np.percentile(d, percentile)
            local_capValue = max(capValue, localPercentile)

            d = np.minimum(d, local_capValue)

            # Trimmed stats
            tMean = np.mean(d).item()
            tErr = np.std(d).item() / np.sqrt(len(d))
            result.append([ref_Name, strand, pos, coverage, tMean, tErr])
    
    print ("one strand done", time.time()-t0)

    for pos in range(start, end):   
        if pos in s1dict:
            d = np.array(s1dict[pos])
            if len(d) <= 2:
                continue

            # NOTE -- this is where the strand flipping occurs -- make sure to
            # reproduce this in the all calling methods
            strand = 1 - 1
            coverage = len(d)

            percentile = min(90, (1.0 - 1.0 / (len(d) - 1)) * 100)
            localPercentile = np.percentile(d, percentile)
            local_capValue = max(capValue, localPercentile)
            d = np.minimum(d, local_capValue)

            # Trimmed stats
            tMean = np.mean(d).item()
            tErr = np.std(d).item() / np.sqrt(len(d))

            result.append([ref_Name, strand, pos, coverage, tMean, tErr])

    print ("raw ipd is counted", time.time()-t0)
    return result

def extract_context(fasta):
    ## load the fasta using biopython
    # print ("loading fasta")
    seq_dict = {}
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        # print(record.id)
        # if record.id != "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L_219069_438138":
        #     continue
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)
        ## convert to capital
        raw_seq = seq.upper()
        seq_dict[record.id] = raw_seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
    return seq_dict

def load_IPD(each_ref, contig_bam, df_file):
    t0 = time.time()
    alignments = AlignmentSet(contig_bam, referenceFastaFname=fasta)
    refInfo = alignments.referenceInfoTable
    # print ("refInfo", refInfo.shape, len(refInfo))
    for my_ref in refInfo:
        if my_ref.Name == each_ref:
            each_ref = my_ref
            break
    print ("ref loaded", each_ref.Name, each_ref.Length, time.time()-t0)
    factor = 1.0 / alignments.readGroupTable[0].FrameRate
    # global alignments
    refGroupId = alignments.referenceInfo(each_ref.Name).Name

    chunks_num = int(each_ref.Length/max_region)
    chunks_num = 1 if chunks_num < 1 else chunks_num
        
    print ("chunks_num", chunks_num)
    result = []
    for i in range(chunks_num):
        start = i*max_region
        end = (i+1)*max_region
        if end > each_ref.Length:
            end = each_ref.Length
        if i == chunks_num-1:
            end = each_ref.Length
        print ("chunk", i, start, end)
        result += _loadRawIpds(alignments, refGroupId, each_ref, start, end, factor, )

    # chunk_result = _loadRawIpds(alignments, refGroupId, each_ref, factor, )
    # print ("rawIpds", rawIpds.shape, time.time()-t0)
    # combined_df = pd.concat(result, ignore_index=True)
    combined_df = pd.DataFrame(result, columns=['refName', 'strand', 'tpl', 'coverage', 'tMean', 'tErr'])
    combined_df.to_csv(df_file, index=False)
    print ("raw ipd df saved", df_file, time.time()-t0)




if __name__ == "__main__":

    # subread_bam = "/home/shuaiw/methylation/data/borg/b_contigs/11.align.bam"
    # fasta = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa"
    # subread_bam = "/home/shuaiw/methylation/data/borg/new_test4/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_267_C.bam"
    # fasta = "/home/shuaiw/methylation/data/borg/new_test4/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_267_C.fa"
    # outputfile = "/home/shuaiw/methylation/data/borg/b_contigs/test/test.csv"

    subread_bam = sys.argv[1]
    fasta = sys.argv[2]
    outputfile = sys.argv[3]

    # subread_bam = "/home/shuaiw/methylation/data/borg/b_contigs/11.align.bam"
    # fasta = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa"
    # outputfile = "/home/shuaiw/methylation/data/borg/b_contigs/test/test_7.csv"

    # subread_bam = "/home/shuaiw/methylation/data/borg/new_test5/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.bam"
    # fasta = "/home/shuaiw/methylation/data/borg/new_test5/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.fa"
    # outputfile = "/home/shuaiw/methylation/data/borg/new_test5/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.ipd_test.csv"


    ## build outdir if not exists
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)

    seq_dict = extract_context(fasta)
    print ("fasta loaded, contig num:", len(seq_dict))

    for each_ref in seq_dict:
        load_IPD(each_ref, subread_bam, outputfile)

    print ("IPD loaded")

    # python /home/shuaiw/Methy/standard_load6.py /home/shuaiw/methylation/data/borg/new_test5/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.bam /home/shuaiw/methylation/data/borg/new_test5/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.fa /home/shuaiw/methylation/data/borg/new_test5/ipd/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_766_C.ipd1.csv
