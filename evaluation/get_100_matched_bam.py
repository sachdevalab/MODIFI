import pysam
import numpy as np
import os


def filter(rawbam, filterbam, max_NM=0): 
    samfile = pysam.AlignmentFile(rawbam, "rb")
    contig_samfile = pysam.AlignmentFile(filterbam, "wb", header=samfile.header)
    valid_num = 0
    for read in samfile:
        ### calculate the alignment identity
        if read.is_unmapped:
            continue
        if not read.has_tag('NM'):
            print (f"Read {read.query_name} has no NM tag")
            continue
        if read.get_tag("NM") > max_NM:
            continue
        ## get cigar, check if it contain insertion or deletion
        cigar = read.cigarstring
        if 'I' in cigar or 'D' in cigar:
            continue

        contig_samfile.write(read)
        valid_num += 1
    print (f"Finished {rawbam}, {valid_num} reads are written")
    samfile.close()
    contig_samfile.close()
    ## index the bam file
    os.system(f"samtools index {filterbam}")  # Index the output BAM file
    return

raw_bam = "/home/shuaiw/methylation/data/borg/bench/break5/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_936_C_1270095_1693460.bam"
filter_bam = "/home/shuaiw/methylation/data/borg/bench/break5/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_936_C_1270095_1693460.filter.bam"
filter(raw_bam, filter_bam)