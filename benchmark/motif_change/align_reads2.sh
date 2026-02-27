#!/bin/bash
#SBATCH --job-name=inversion     # Job name
#SBATCH --partition=standard # Partition name

## construct output directory if not exists


outdir=/home/shuaiw/borg/paper/motif_change/new_alignment2/
ref=/home/shuaiw/borg/paper/run2/infant_14/infant_14_methylation4/contigs/infant_14_31_C.fa

if [ ! -d $outdir ]; then
    mkdir $outdir
fi


ccs_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_2G1_pacbio.bam
prefix=$outdir/week2
align_bam=$prefix.align.ccs.bam



/usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref $ccs_bam $align_bam \
 --preset CCS -j $SLURM_CPUS_ON_NODE --sort -J $SLURM_CPUS_ON_NODE
samtools index $align_bam
/home/shuaiw//smrtlink/pbindex $align_bam
samtools view -b $align_bam infant_14_31_C:2273856-2377115 > $prefix.region.bam

ccs_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_3G1_pacbio.bam
prefix=$outdir/week3
align_bam=$prefix.align.ccs.bam


/usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref $ccs_bam $align_bam \
 --preset CCS -j $SLURM_CPUS_ON_NODE --sort -J $SLURM_CPUS_ON_NODE
samtools index $align_bam
/home/shuaiw//smrtlink/pbindex $align_bam
samtools view -b $align_bam infant_14_31_C:2273856-2377115 > $prefix.region.bam


ccs_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_4G1_pacbio.bam
prefix=$outdir/week4
align_bam=$prefix.align.ccs.bam


/usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref $ccs_bam $align_bam \
 --preset CCS -j $SLURM_CPUS_ON_NODE --sort -J $SLURM_CPUS_ON_NODE
samtools index $align_bam
/home/shuaiw//smrtlink/pbindex $align_bam
samtools view -b $align_bam infant_14_31_C:2273856-2377115 > $prefix.region.bam



ccs_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_5G1_pacbio.bam
prefix=$outdir/week5
align_bam=$prefix.align.ccs.bam


/usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref $ccs_bam $align_bam \
 --preset CCS -j $SLURM_CPUS_ON_NODE --sort -J $SLURM_CPUS_ON_NODE
samtools index $align_bam
/home/shuaiw//smrtlink/pbindex $align_bam
samtools view -b $align_bam infant_14_31_C:2273856-2377115 > $prefix.region.bam



ccs_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_6G1_pacbio.bam
prefix=$outdir/week6
align_bam=$prefix.align.ccs.bam


/usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref $ccs_bam $align_bam \
 --preset CCS -j $SLURM_CPUS_ON_NODE --sort -J $SLURM_CPUS_ON_NODE
samtools index $align_bam
/home/shuaiw//smrtlink/pbindex $align_bam
samtools view -b $align_bam infant_14_31_C:2273856-2377115 > $prefix.region.bam
