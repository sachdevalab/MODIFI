import pysam
import random

def downsample_subreads_bam(input_bam, output_bam, fraction=0.01, seed=42):
    """
    Randomly downsample a PacBio subreads BAM file.
    """
    random.seed(seed)
    
    with pysam.AlignmentFile(input_bam, "rb", check_sq=False) as in_bam:
        with pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
            for read in in_bam:
                if random.random() < fraction:
                    out_bam.write(read)
# Example usage
downsample_subreads_bam(
    input_bam="/groups/banfield/projects/environmental/sr/srvp2022/pacbio/raw.d/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI/download/XRSBK_20221007_S64018_PL100268287-1_C01.subreads.bam",
    output_bam="/home/shuaiw/borg/allison/ecoli/subsampled.subreads.bam",
    fraction=0.01
)

