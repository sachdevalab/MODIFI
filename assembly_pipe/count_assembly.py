"""
count how many reads can map to assembly
"""


import os
import gzip
import pysam
import subprocess
import sys
import pandas as pd


def count_reads(align_bam, ccs_bam, align_count):
    """
    ## count number of reads in raw fastq in python
    """

    # Count the number of reads in the raw ccs bam using samtools
    

    result = subprocess.run(
        ["samtools", "view", "-c", ccs_bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )
    raw_read_count = int(result.stdout.strip())
    print(f"Total reads: {result.stdout.strip()}")

    # Open the BAM file
    bamfile = pysam.AlignmentFile(align_bam, "rb")

    # Count the number of reads in the BAM file
    # read_count = sum(1 for read in bamfile.fetch(until_eof=True) if not read.is_unmapped)
    read_set = set()
    for read in bamfile.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        ## map quality should be greater than 20
        if read.mapping_quality >= 20:
            read_set.add(read.query_name)
    read_count = len(read_set)
    bamfile.close()




    ## count the ratio of reads that can map to assembly
    if raw_read_count == 0:
        ratio = 0
    else:
        ratio = read_count / raw_read_count
    print (f"Read count: {read_count}, Raw read count: {raw_read_count}, Ratio: {ratio:.4f}")
    data = {
        "map_read_count": read_count,
        "raw_read_count": raw_read_count,
        "ratio": ratio
    }
    ## convert to df and output
    
    df = pd.DataFrame([data])
    df.to_csv(align_count, index=False)

    return read_count, raw_read_count, ratio


if __name__ == "__main__":
    # align_bam = "/home/shuaiw/borg/paper/run/ERR12723529/ERR12723529.align.bam"
    # ccs_bam = "/home/shuaiw/borg/paper/aws/ERR12723528/ERR12723528.ccs.bam"

    prefix = sys.argv[1] 
    work_dir = sys.argv[2]
    ccs_bam = sys.argv[3]

    align_bam = os.path.join(work_dir, f"{prefix}.align.bam")
    align_count = os.path.join(work_dir, f"{prefix}.align.count.csv")
    count_reads(align_bam, ccs_bam, align_count)
