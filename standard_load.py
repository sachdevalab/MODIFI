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

# subread_bam = "/home/shuaiw/methylation/data/borg/human/human_000733.subreads.align.bam"

# subread_bam = "/home/shuaiw/methylation/data/borg/split_bam_dir2/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_7580_L.bam"
# fasta = "/home/shuaiw/methylation/data/borg/split_bam_dir2/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_7580_L.fasta"

subread_bam = "/home/shuaiw/methylation/data/borg/b_contigs/1.align.bam"
fasta = "/home/shuaiw/methylation/data/borg/b_contigs/contigs/1.fa"

print ("loading alignments", subread_bam)
alignments = AlignmentSet(subread_bam,
                                referenceFastaFname=fasta)

refInfo = alignments.referenceInfoTable
print ("refInfo", refInfo)


refGroupId = alignments.referenceInfo(
        'SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C_0_852595').Name
hits = [hit for hit in alignments.readsInRange(refGroupId,
                                                    max(0, 0), 10000)]

rawIpds = _loadRawIpds(hits, 0, 10000, 1.0)
print ("rawIpds", rawIpds)
for i in range(10):
    print (rawIpds[i])

ipdVect = rawIpds['ipd']
# print (ipdVect.mean().item(), rawIpds.size)
## cal the mean of the locus 1 and strand 0
first_locus = []
for i in range(len(rawIpds)):
    if rawIpds[i]['tpl'] == 2 and rawIpds[i]['strand'] == 0:
        print (rawIpds[i])
        # print ("")
        first_locus.append(rawIpds[i]['ipd'])
print ("mean of the locus 1 and strand 0", np.mean(first_locus), len(first_locus)  )


if ipdVect.size < 10:
    # Default is there is no coverage
    capValue = 5.0
else:
    # Compute IPD quantiles on the current block -- will be used for
    # trimming extreme IPDs
    capValue = np.percentile(ipdVect, 99)

print ("capValue", capValue)