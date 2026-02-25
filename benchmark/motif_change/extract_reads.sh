bam=/home/shuaiw/borg/paper/run2/infant_2/infant_2_methylation4/bams/infant_2_3_C.bam
samtools view -b $bam infant_2_3_C:2889998-2893008 > region_reads.bam