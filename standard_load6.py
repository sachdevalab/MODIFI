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

# alignments = None

# Raw ipd record
ipdRec = [('tpl', '<u4'), ('strand', '<i8'), ('ipd', '<f4')]

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

def _loadRawIpds(alignments, refGroupId, each_ref, df_file, factor=1.0):
    """
    Get a DataFrame of the raw ipds in the give alignment hits, indexed by template position and strand.
    Factor is a normalization factor to the get units into seconds.
    """

    # Put in an empty 'starter' array -- the np.concatenate call below will
    # fail on an empty list
    array0 = np.zeros(0, dtype=ipdRec)
    t0 = time.time()

    # Maintain separate lists for each strand to speed up sorting
    s0dict = defaultdict(list)
    s1dict = defaultdict(list)
    ipdVect = []

    for aln in alignments.readsInRange(refGroupId, 0, each_ref.Length):
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
        np.logical_and(referencePositions < each_ref.Length, matched, matched)
        np.logical_and(referencePositions >= 0, matched, matched)
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

        # dfTemp = np.zeros(nm, dtype=ipdRec)
        # dfTemp['ipd'] = ipd
        # dfTemp['tpl'] = tpl
        # dfTemp['strand'] = aln.isReverseStrand
        # print ("ipd", ipd)
        # ipdVect = np.concatenate([ipdVect, ipd])
        # if aln.isForwardStrand:
        #     for i in range(ipd.size):
        #         s0dict[tpl[i]].append(ipd[i])
        # else:
        #     for i in range(ipd.size):
        #         s1dict[tpl[i]].append(ipd[i])
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

    # s0Ipds, s1Ipds = [], []
    result = []
    for pos in range(0, each_ref.Length):
        if pos in s0dict:
            res = dict()
            d = np.array(s0dict[pos])
            if d.size <= 2:
                continue

            res['refName'] = each_ref.Name

            # NOTE -- this is where the strand flipping occurs -- make sure to
            # reproduce this in the all calling methods
            strand = res['strand'] = 1 - 0
            res['tpl'] = pos
            res['coverage'] = d.size

            percentile = min(90, (1.0 - 1.0 / (d.size - 1)) * 100)
            localPercentile = np.percentile(d, percentile)
            local_capValue = max(capValue, localPercentile)

            # np.minimum(d, capValue, out=d)  # this version will send capped IPDs
            # to modified fraction estimator
            d = np.minimum(d, local_capValue)

            # Trimmed stats
            res['tMean'] = d.mean().item()
            res['tErr'] = np.std(d).item() / np.sqrt(d.size)

            result.append(res)
    
    print ("one strand done", time.time()-t0)

    for pos in range(0, each_ref.Length):   
        if pos in s1dict:
            res = dict()
            d = np.array(s1dict[pos])
            if d.size <= 2:
                continue

            res['refName'] = each_ref.Name

            # NOTE -- this is where the strand flipping occurs -- make sure to
            # reproduce this in the all calling methods
            strand = res['strand'] = 1 - 1
            res['tpl'] = pos
            res['coverage'] = d.size

            percentile = min(90, (1.0 - 1.0 / (d.size - 1)) * 100)
            localPercentile = np.percentile(d, percentile)
            local_capValue = max(capValue, localPercentile)

            # np.minimum(d, capValue, out=d)  # this version will send capped IPDs
            # to modified fraction estimator
            d = np.minimum(d, local_capValue)

            # Trimmed stats
            res['tMean'] = d.mean().item()
            res['tErr'] = np.std(d).item() / np.sqrt(d.size)

            result.append(res)

    print ("raw ipd is counted", time.time()-t0)
    # combined_df = pd.concat(result, ignore_index=True)
    combined_df = pd.DataFrame(result)
    combined_df.to_csv(df_file, index=False)
    print ("raw ipd df saved", df_file)

    # for pos in range(0, each_ref.Length):
    #     if pos in s1dict:
    #         obj = {'tpl': pos, 'strand': "1",
    #                 'data': np.array(s1dict[pos])}
    #         views.append(obj)
        
    # print ("start sorting")
    # print (s0list)
    # t0 = time.time()
    # # Sort the set of ipd observations
    # s0Ipds = np.concatenate(s0list)
    # sortOrder = np.argsort(s0Ipds['tpl'])
    # s0Ipds = s0Ipds[sortOrder]

    # s1Ipds = np.concatenate(s1list)
    # sortOrder = np.argsort(s1Ipds['tpl'])
    # s1Ipds = s1Ipds[sortOrder]
    # print ("sorting done", time.time()-t0)
    # print (s1Ipds)


def _loadRawIpds2(alignments, refGroupId, each_ref, factor=1.0):
    """
    Get a DataFrame of the raw ipds in the give alignment hits, indexed by template position and strand.
    Factor is a normalization factor to the get units into seconds.
    """

    # Put in an empty 'starter' array -- the np.concatenate call below will
    # fail on an empty list
    array0 = np.zeros(0, dtype=ipdRec)

    # Maintain separate lists for each strand to speed up sorting
    s0list = [array0]
    s1list = [array0]

    for aln in alignments.readsInRange(refGroupId, 0, each_ref.Length):
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
        np.logical_and(referencePositions < each_ref.Length, matched, matched)
        np.logical_and(referencePositions >= 0, matched, matched)
        nm = matched.sum()

        # Bail out if we don't have any samples
        if nm == 0:
            continue

        ipd = rawIpd[matched]
        tpl = referencePositions[matched]

        dfTemp = np.zeros(nm, dtype=ipdRec)
        dfTemp['ipd'] = ipd
        dfTemp['tpl'] = tpl
        dfTemp['strand'] = aln.isReverseStrand

        if aln.isForwardStrand:
            s0list.append(dfTemp)
        else:
            s1list.append(dfTemp)

    print ("start sorting")
    print (s0list)
    t0 = time.time()
    # Sort the set of ipd observations
    s0Ipds = np.concatenate(s0list)
    sortOrder = np.argsort(s0Ipds['tpl'])
    s0Ipds = s0Ipds[sortOrder]

    s1Ipds = np.concatenate(s1list)
    sortOrder = np.argsort(s1Ipds['tpl'])
    s1Ipds = s1Ipds[sortOrder]
    print ("sorting done", time.time()-t0)
    print (s1Ipds)
    return np.concatenate([s0Ipds, s1Ipds])

def _chunkRawIpds(rawIpds):
    """
    Return a list of view recarrays into the rawIpds recarray, one for each unique (tpl, stand) level
    """
    views = []

    # Bail out if we have no data
    if rawIpds.size == 0:
        return views

    start = 0
    tpl = rawIpds['tpl']
    strand = rawIpds['strand']

    # Start off at the first chunk
    curIdx = (tpl[0], strand[0])
    for i in range(1, rawIpds.shape[0]):
        newIdx = (tpl[i], strand[i])

        # In this case we are still int he same chunk -- continue
        if curIdx == newIdx:
            continue

        # In this case we have completed the chunk -- emit the chunk
        else:
            obj = {'tpl': curIdx[0], 'strand': curIdx[1],
                    'data': rawIpds[start:i]}
            views.append(obj)
            start = i
            curIdx = newIdx
    # Make sure to return final chunk
    obj = {'tpl': curIdx[0], 'strand': curIdx[1], 'data': rawIpds[start:]}
    views.append(obj)

    # # If the user has specified a maximum coverage level to use, enforce it
    # # here -- just take the first n reads
    # if self.options.maxCoverage is not None:
    #     maxCov = self.options.maxCoverage
    #     for x in views:
    #         d = x['data']
    #         d = d[0:maxCov]
    #         x['data'] = d

    return views

def _computePositionSyntheticControl(caseObservations, capValue, refId, refName):
    """Summarize the observed ipds at one template position/strand, using the synthetic ipd model"""

    # Compute stats on the observed ipds
    d = caseObservations['data']['ipd']
    res = dict()

    # ref00000x name
    # res['refId'] = refId

    # FASTA header name
    res['refName'] = refName

    # NOTE -- this is where the strand flipping occurs -- make sure to
    # reproduce this in the all calling methods
    strand = res['strand'] = 1 - caseObservations['strand']
    tpl = res['tpl'] = caseObservations['tpl']
    res['coverage'] = d.size

    # Don't compute these stats - they just take time and confuse things
    # res['mean'] = d.mean().item()
    # res['median'] = np.median(d).item()
    # res['std'] = np.std(d).item()
    # Compute the predicted IPD from the model
    # NOTE! The ipd model is in the observed read strand
    # if modelPrediction is None:
    #     modelPrediction = self.meanIpdFunc(tpl, strand).item()
    # res['modelPrediction'] = modelPrediction

    # res['base'] = self.cognateBaseFunc(tpl, strand)

    # Store in case of methylated fraction estimtion:
    # res['rawData'] = d   # discard it

    percentile = min(90, (1.0 - 1.0 / (d.size - 1)) * 100)
    localPercentile = np.percentile(d, percentile)
    capValue = max(capValue, localPercentile)

    # np.minimum(d, capValue, out=d)  # this version will send capped IPDs
    # to modified fraction estimator
    d = np.minimum(d, capValue)

    # Trimmed stats
    res['tMean'] = d.mean().item()
    res['tErr'] = np.std(d).item() / np.sqrt(d.size)
    # print (res['tpl'], res['strand'], res['tMean'], res['tErr'], res['coverage'])


    return res
      
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
    # hits = [hit for hit in alignments.readsInRange(refGroupId, 0, each_ref.Length)]
    # hits = [hit for hit in alignments.readsInRange(refGroupId, 0, 1000)]
    print ("hits", time.time()-t0)

    _loadRawIpds(alignments, refGroupId, each_ref, df_file, factor, )
    # print ("rawIpds", rawIpds.shape, time.time()-t0)
    """
    caseChunks = _chunkRawIpds(rawIpds)
    print ("chunks", time.time()-t0)

    ipdVect = rawIpds['ipd']
    if ipdVect.size < 10:
        # Default is there is no coverage
        capValue = 5.0
    else:
        # Compute IPD quantiles on the current block -- will be used for
        # trimming extreme IPDs
        capValue = np.percentile(ipdVect, 99)
    print ("capValue", capValue, time.time()-t0)

    # goodSites = [x for x in caseChunks if x['data']['ipd'].size > 2]
    result = []
    for x in caseChunks:
        if x['data']['ipd'].size > 2:
            pass
        else:
            continue
        # print (x['tpl'], x['strand'], x['data']['ipd'].size)
        res = _computePositionSyntheticControl(x, capValue, refGroupId, each_ref.Name)
        # if res['strand'] == 1:
        result.append(pd.DataFrame([res]))

    print ("raw ipd is counted", time.time()-t0)
    combined_df = pd.concat(result, ignore_index=True)
    # print ("combined_df", combined_df)
    print (len(combined_df), 'rows')

    # df_file = os.path.join(outdir, each_ref.Name + ".ipd1.csv")
    combined_df.to_csv(df_file, index=False)
    print ("raw ipd df saved", df_file)
    return combined_df
    """



if __name__ == "__main__":
    up = 8
    down = 4

    subread_bam = "/home/shuaiw/methylation/data/borg/b_contigs/11.align.bam"
    fasta = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa"
    outputfile = "/home/shuaiw/methylation/data/borg/b_contigs/test/test.csv"

    # subread_bam = sys.argv[1]
    # fasta = sys.argv[2]
    # outputfile = sys.argv[3]

    # subread_bam = "/home/shuaiw/borg/break_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
    # fasta = "/home/shuaiw/borg/contigs/break_contigs.fasta"
    # outdir = "/home/shuaiw/methylation/data/borg/b_contigs/ipds8/"
    # bam_dir = "/home/shuaiw/methylation/data/borg/b_contigs/bams/"

    # subread_bam = "/home/shuaiw/borg/large_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
    # fasta = "/home/shuaiw/borg/contigs/large.contigs.fa"
    # outdir = "/home/shuaiw/methylation/data/borg/large_contigs/ipds/"


    ## build outdir if not exists
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)

    seq_dict = extract_context(fasta)
    print ("fasta loaded, contig num:", len(seq_dict))

    for each_ref in seq_dict:
        load_IPD(each_ref, subread_bam, outputfile)

    print ("IPD loaded")
