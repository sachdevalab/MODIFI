"""
count how many reads can map to assembly
"""


import os


def count_reads(align_bam, raw_fq):
    """
    ## count number of reads in raw fastq in python
    """
    import pysam

    # Open the BAM file
    bamfile = pysam.AlignmentFile(align_bam, "rb")

    # Count the number of reads in the BAM file
    read_count = bamfile.count(until_eof=True, read_callback='all', mapped=True)

    # Close the BAM file
    bamfile.close()

    # Count the number of reads in the raw fastq file
    with open(raw_fq, 'r') as fq:
        raw_read_count = sum(1 for line in fq if line.startswith('@')) // 4

    ## count the ratio of reads that can map to assembly
    if raw_read_count == 0:
        ratio = 0
    else:
        ratio = read_count / raw_read_count
    print (f"Read count: {read_count}, Raw read count: {raw_read_count}, Ratio: {ratio:.4f}")
    return read_count, raw_read_count, ratio


if __name__ == "__main__":
    align_bam = "/home/shuaiw/borg/paper/run/ERR12723529/ERR12723529.align.bam"
    raw_fq = "/home/shuaiw/borg/paper/run/ERR12723529/ERR12723529.hifi.qc.HR.fq.gz"
    count_reads(align_bam, raw_fq)