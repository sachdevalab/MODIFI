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
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager, Process, Lock



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

def _loadRawIpds(alnHitIter, targetStart=-
                     1, targetEnd=3e12, factor=1.0):
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

    for aln in alnHitIter:
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
        np.logical_and(referencePositions < targetEnd, matched, matched)
        np.logical_and(referencePositions >= targetStart, matched, matched)
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

    # Sort the set of ipd observations
    s0Ipds = np.concatenate(s0list)
    sortOrder = np.argsort(s0Ipds['tpl'])
    s0Ipds = s0Ipds[sortOrder]

    s1Ipds = np.concatenate(s1list)
    sortOrder = np.argsort(s1Ipds['tpl'])
    s1Ipds = s1Ipds[sortOrder]

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

def get_IPD_list(contig_length, contig_dict):
    average_IPD = []
    for pos in range(contig_length):
        if pos not in contig_dict:
            average_IPD.append(0)
        else:
            average_IPD.append(round(np.mean(contig_dict[pos]), 1))
    return average_IPD
             
def extract_context(fasta):
    ## load the fasta using biopython
    print ("loading fasta")
    seq_dict = {}
    from Bio import SeqIO
    for record in SeqIO.parse(fasta, "fasta"):
        print(record.id)
        # if record.id != "NC_000001.11":
        #     continue
        ## convert the sequence to string of number, 0 for A, 1 for C, 2 for G, 3 for T, 4 for N
        seq = str(record.seq)
        ## convert to capital
        raw_seq = seq.upper()
        seq_dict[record.id] = raw_seq
        # seq = raw_seq.replace('A', '0').replace('C', '1').replace('G', '2').replace('T', '3').replace('N', '4')
    return seq_dict

def prepare_data(seq, ipd_list, kmer_baseline_dict, kmer_num_dict, up=7, down=3):
    # save_kmer = {}  ## {pos:kmer}
    
    # y = control_list[up:len(seq) - down]
    for i in range(up, len(seq) - down):
        kmer = seq[i-up:i+down]
        if 'N' in kmer:
            continue  # Skip kmers containing 'N'
        ipd = float(ipd_list[i])
        if ipd == 0:
            continue
        # kmer_baseline_dict[kmer].append(ipd)
        kmer_baseline_dict[kmer] += ipd
        kmer_num_dict[kmer] += 1
        # save_kmer[i] = kmer    
    return kmer_baseline_dict, kmer_num_dict

class Contig:
    def __init__(self, name):
        self.name = name
        self.forward_dict = {}
        self.reverse_dict = {}
        self.forward_ipd_sum = {}
        self.reverse_ipd_sum = {}

def get_reverse_cmplement(kmer):
    kmer_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    complement_kmer = ''
    for i in range(len(kmer)):
        complement_kmer += kmer_dict[kmer[i]]
    # reverse_kmer = complement_kmer[::-1]
    # return reverse_kmer
    return complement_kmer

def count_kmer(contig_forward_dict_dict, seq, strand = 1):

    # kmer_baseline_dict = defaultdict(list)
    kmer_baseline_dict = defaultdict(float)   ## sum
    kmer_num_dict = defaultdict(int)  # count
    for contig in contig_forward_dict_dict:
        
        contig_forward_dict = contig_forward_dict_dict[contig]
        # print (contig, len(contig_forward_dict),contig_forward_dict)

        seq = seq_dict[contig]
        if strand == 0:
            seq = get_reverse_cmplement(seq)
        observed_IPD_list = get_IPD_list(len(seq), contig_forward_dict)
        kmer_baseline_dict, kmer_num_dict = prepare_data(seq, observed_IPD_list, kmer_baseline_dict, kmer_num_dict, up, down)
    print ("kmer is counted")
    
    mean_dict, median_dict = {}, {}
    for kmer in kmer_baseline_dict:
        mean_dict[kmer] = round(kmer_baseline_dict[kmer]/kmer_num_dict[kmer], 3) ## calculate the mean of each kmer
        median_dict[kmer] = 0
    print ("mean and median is computed", len(mean_dict), 'kmers')
    return mean_dict, median_dict, kmer_baseline_dict, kmer_num_dict

def align_kmer(contig_forward_dict_dict, \
               kmer_baseline_dict, kmer_num_dict, mean_dict, median_dict, df, strand = 1):
    
    for contig in contig_forward_dict_dict:
        contig_forward_dict = contig_forward_dict_dict[contig]
        print (contig, len(contig_forward_dict))

        seq = seq_dict[contig]
        for pos in range(up, len(seq) - down):
            kmer = seq[pos-up:pos+down]
            if 'N' in kmer:
                continue  # Skip kmers containing 'N'
            # if kmer_num_dict[kmer] < 10:
            #     continue
            if kmer in kmer_baseline_dict and pos in contig_forward_dict:
                # print (pos, kmer, len(kmer_baseline_dict[kmer]), np.mean(kmer_baseline_dict[kmer]), contig_forward_dict[pos][0], ipd_sum_for_control[pos])

                ## To do, make it more quick
                row_index = df.loc[(df['tpl'] == pos) & (df['strand'] == strand)].index
                if not row_index.empty:
                    df.loc[row_index, 'kmer'] = kmer
                    df.loc[row_index, 'count'] = kmer_num_dict[kmer]
                    df.loc[row_index, 'mean'] = mean_dict[kmer]
                    df.loc[row_index, 'median'] = median_dict[kmer]
    return df

def get_kmer(seq, i, up=7, down=3):
    if i - up < 0 or i + down > len(seq):
        return None
    kmer = seq[i-up:i+down]
    if 'N' in kmer:
        return None
    return kmer

def load_IPD1(each_ref):
    print ('hi')
    # t0 = time.time()
    alignments = AlignmentSet(subread_bam, referenceFastaFname=fasta)
    refGroupId = alignments.referenceInfo(each_ref.Name).Name
    for hit in alignments.readsInRange(refGroupId, 0, each_ref.Length):
        pass
    # hits = [hit for hit in alignments.readsInRange(refGroupId, 0, each_ref.Length)]
    # hits = [hit for hit in alignments.readsInRange(refGroupId, 0, 100000)]
    print ("11 hits", alignments.referenceInfo(each_ref.Name).Name)

def load_IPD(each_ref, kmer_baseline_dict, kmer_num_dict):
    # refGroupId = alignments.referenceInfo(
    #         'SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C_0_852595').Name
    t0 = time.time()
    alignments = AlignmentSet(subread_bam, referenceFastaFname=fasta)
    refGroupId = alignments.referenceInfo(each_ref.Name).Name
    hits = [hit for hit in alignments.readsInRange(refGroupId, 0, each_ref.Length)]
    # hits = [hit for hit in alignments.readsInRange(refGroupId, 0, 1000)]
    print ("hits", len(hits), time.time()-t0)


    rawIpds = _loadRawIpds(hits, 0, each_ref.Length, factor)
    print ("rawIpds", rawIpds.shape)
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
    seq = seq_dict[each_ref.Name]
    complement_seq = get_reverse_cmplement(seq)
    print ("goodSites", time.time()-t0)


    # for x in goodSites:
    for x in caseChunks:
        if x['data']['ipd'].size > 2:
            pass
        else:
            continue
        # print (x['tpl'], x['strand'], x['data']['ipd'].size)
        res = _computePositionSyntheticControl(x, capValue, refGroupId, each_ref.Name)

        
        if res['strand'] == 1:
            kmer = get_kmer(seq, res['tpl'])
            res['kmer'] = kmer
        else:
            res['kmer'] = None
            continue
        
        if kmer is None:
            continue

        # kmer_baseline_dict[kmer] += res['tMean']
        # kmer_num_dict[kmer] += 1
        with lock:
            if kmer not in kmer_baseline_dict:
                kmer_baseline_dict[kmer] = 0.0
            kmer_baseline_dict[kmer] += res['tMean']
            if kmer not in kmer_num_dict:
                kmer_num_dict[kmer] = 0
            kmer_num_dict[kmer] += 1

        result.append(pd.DataFrame([res]))
    print ("kmer is counted", len(kmer_baseline_dict), len(kmer_num_dict), time.time()-t0)

    
    combined_df = pd.concat(result, ignore_index=True)
    # print ("combined_df", combined_df)
    print (len(combined_df), 'rows')

    df_file = os.path.join(outdir, each_ref.Name + ".ipd1.csv")
    combined_df.to_csv(df_file, index=False)
    print ("raw ipd df saved", df_file)

    # all_ref_IPDs[each_ref.Name] = combined_df
    return combined_df

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
    # print (x, mu, sigma, p_value)
    return p_value

def normalize_IPD(each_ref, combined_df, kmer_mean_dict, kmer_num_dict):
    combined_df['count'] = None
    combined_df['control'] = None

    for index, row in combined_df.iterrows():
        kmer = row['kmer']
        if kmer is None:
            combined_df.loc[index, 'count'] = 0
            combined_df.loc[index, 'control'] = 1
        else: 
            combined_df.loc[index, 'count'] = kmer_num_dict[kmer]
            combined_df.loc[index, 'control'] = kmer_mean_dict[kmer]

    combined_df = get_ipd_ratio(combined_df)
    df_file = os.path.join(outdir, each_ref.Name + ".ipd2.csv")
    combined_df.to_csv(df_file, index=False)
    print ("df saved", df_file)

def get_ipd_ratio(df):
    print ('get ipd_ratio')
    ## get ratio which is the ratio between estimated and ipdsum
    df['ipd_ratio'] = df['tMean'] / df['control'] 
    # data = df['ipd_ratio'].values.reshape(-1, 1)
    # gmm = GaussianMixture(n_components=1, covariance_type='full', random_state=42)
    # gmm.fit(data)
    # mean = gmm.means_[0][0]
    # std = np.sqrt(gmm.covariances_[0][0])
    ## simply calculate mean and std for ipd_ratio
    mean = df['ipd_ratio'].mean()
    std = df['ipd_ratio'].std()
    ## ofr each value, calculate the probability belong to the model
    ## add pvalue column to the dataframe
    df['pvalue'] = df['ipd_ratio'].apply(lambda x: p_value_right_tail(x, mean, std))
    return df

def process_reference(each_ref):
    print(f"{each_ref.Name}, {each_ref.Length}")
    return load_IPD(each_ref)

def normalize_IPD_wrapper(ref):
    try:
        df_file = os.path.join(outdir, ref.Name + ".ipd1.csv")
        combined_df = pd.read_csv(df_file)
        normalize_IPD(ref, combined_df, kmer_mean_dict, kmer_num_dict)
        return ref.Name
    except Exception as e:
        print(f"An error occurred while processing {ref.Name}: {e}")
        return None

if __name__ == "__main__":
    up = 8
    down = 4

    subread_bam = "/home/shuaiw/methylation/data/borg/b_contigs/12.align.bam"
    fasta = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/12.fa"
    outdir = "/home/shuaiw/methylation/data/borg/b_contigs/"


    # subread_bam = "/home/shuaiw/borg/break_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
    # fasta = "/home/shuaiw/borg/contigs/break_contigs.fasta"
    # outdir = "/home/shuaiw/methylation/data/borg/b_contigs/ipds2/"

    num_processes = 10

    ## build outdir if not exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print ("loading alignments", subread_bam)
    alignments = AlignmentSet(subread_bam,
                                    referenceFastaFname=fasta)
    factor = 1.0 / alignments.readGroupTable[0].FrameRate  ## The frame rate represents the speed of data acquisition in frames per second during sequencing.
    print ("factor", factor)

    refInfo = alignments.referenceInfoTable
    print ("refInfo", refInfo.shape, len(refInfo), refInfo)

    seq_dict = extract_context(fasta)
    print ("ref loaded", len(seq_dict))

    manager = Manager()
    kmer_baseline_dict = manager.dict()
    kmer_num_dict = manager.dict()
    lock = Lock()

    # kmer_baseline_dict = defaultdict(float)   ## sum
    # kmer_num_dict = defaultdict(int)  # count

    all_ref_IPDs = {}

    # for each_ref in refInfo:
    #     print (each_ref.Name, each_ref.Length)
    #     combined_df = load_IPD(each_ref, kmer_baseline_dict, kmer_num_dict)

    # Using ProcessPoolExecutor for multithreading
    with ProcessPoolExecutor(max_workers=num_processes) as executor:  # Adjust max_workers as needed
        future_to_ref = {executor.submit(load_IPD, ref, kmer_baseline_dict, kmer_num_dict): ref for ref in refInfo}
        results = []
        for future in as_completed(future_to_ref):
            ref = future_to_ref[future]
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"An error occurred while processing {ref}: {e}")

    
    print ("IPD loaded")
    print ("kmer saved", len(kmer_baseline_dict), len(kmer_num_dict))
    ## cal mean for each kmer
    kmer_mean_dict = {}
    for kmer in kmer_baseline_dict:
        kmer_mean_dict[kmer] = round(kmer_baseline_dict[kmer]/kmer_num_dict[kmer], 3)
    print ("mean is computed", len(kmer_mean_dict), 'kmers')
    ## use pickle to save the kmer_mean_dict
    with open(os.path.join(outdir, 'kmer_mean_dict.pkl'), 'wb') as f:
        pickle.dump(kmer_mean_dict, f)
    ## also save kmer_num_dict
    with open(os.path.join(outdir, 'kmer_num_dict.pkl'), 'wb') as f:
        pickle.dump(kmer_num_dict, f)

    ## load kmer_mean_dict
    # with open(os.path.join(outdir, 'kmer_mean_dict.pkl'), 'rb') as f:
    #     kmer_mean_dict = pickle.load(f)
    #     print ("kmer_mean_dict loaded", len(kmer_mean_dict))

    # for each_ref in refInfo:
    #     # combined_df = all_ref_IPDs[each_ref.Name]
    #     df_file = os.path.join(outdir, each_ref.Name + ".ipd1.csv")
    #     combined_df = pd.read_csv(df_file)
    #     normalize_IPD(each_ref, combined_df, kmer_mean_dict, kmer_num_dict)

    # Using ProcessPoolExecutor for multithreading
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_ref = {executor.submit(normalize_IPD_wrapper, ref): ref for ref in refInfo}
        for future in as_completed(future_to_ref):
            ref = future_to_ref[future]
            try:
                result = future.result()
                if result is not None:
                    print(f"Normalization completed for {result}")
            except Exception as e:
                print(f"An error occurred while processing {ref.Name}: {e}")




