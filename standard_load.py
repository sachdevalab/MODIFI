"""
load ipd from bam, refer to
https://github.com/PacificBiosciences/kineticsTools/blob/master/kineticsTools/KineticWorker.py#L661
"""


from pbcore.io import AlignmentSet
import numpy as np


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

def _computePositionSyntheticControl(caseObservations, capValue):
    """Summarize the observed ipds at one template position/strand, using the synthetic ipd model"""

    # Compute stats on the observed ipds
    d = caseObservations['data']['ipd']
    res = dict()

    # # ref00000x name
    # res['refId'] = self.refId

    # # FASTA header name
    # res['refName'] = self.refName

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
    res['rawData'] = d

    percentile = min(90, (1.0 - 1.0 / (d.size - 1)) * 100)
    localPercentile = np.percentile(d, percentile)
    capValue = max(capValue, localPercentile)

    # np.minimum(d, capValue, out=d)  # this version will send capped IPDs
    # to modified fraction estimator
    d = np.minimum(d, capValue)

    # Trimmed stats
    res['tMean'] = d.mean().item()
    res['tErr'] = np.std(d).item() / np.sqrt(d.size)
    print (res['tpl'], res['strand'], res['tMean'], res['tErr'], res['coverage'])

# subread_bam = "/home/shuaiw/methylation/data/borg/human/human_000733.subreads.align.bam"

# subread_bam = "/home/shuaiw/methylation/data/borg/split_bam_dir2/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_7580_L.bam"
# fasta = "/home/shuaiw/methylation/data/borg/split_bam_dir2/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_7580_L.fasta"

subread_bam = "/home/shuaiw/methylation/data/borg/b_contigs/1.align.bam"
fasta = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/1.fa"

print ("loading alignments", subread_bam)
alignments = AlignmentSet(subread_bam,
                                referenceFastaFname=fasta)

refInfo = alignments.referenceInfoTable
print ("refInfo", refInfo.shape, len(refInfo), refInfo)


refGroupId = alignments.referenceInfo(
        'SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C_0_852595').Name
hits = [hit for hit in alignments.readsInRange(refGroupId, 0, 10000)]

factor = 1.0 / alignments.readGroupTable[0].FrameRate  ## The frame rate represents the speed of data acquisition in frames per second during sequencing.
print ("factor", factor)
rawIpds = _loadRawIpds(hits, 0, 1000, factor)
print ("rawIpds", rawIpds.shape)
caseChunks = _chunkRawIpds(rawIpds)
# print ("chunks", chunks)

ipdVect = rawIpds['ipd']
if ipdVect.size < 10:
    # Default is there is no coverage
    capValue = 5.0
else:
    # Compute IPD quantiles on the current block -- will be used for
    # trimming extreme IPDs
    capValue = np.percentile(ipdVect, 99)

print ("capValue", capValue)

goodSites = [x for x in caseChunks if x['data']['ipd'].size > 2]
for x in goodSites:
    print (x['tpl'], x['strand'], x['data']['ipd'].size)
    _computePositionSyntheticControl(x, capValue)
    # break
    if x['tpl'] > 5:
        break


