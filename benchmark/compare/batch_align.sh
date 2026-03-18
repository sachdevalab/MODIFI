#!/bin/bash
#SBATCH --job-name=ipd_28      # Job name
#SBATCH --partition=standard # Partition name

subreads_bam=/home/shuaiw/borg/XRSBK_20221007_S64018_PL100268287-1_C01.subreads.bam

ref_dir=/home/shuaiw/borg/contigs/

ref_set=(test_100.fa test_200.fa test_300.fa test_500.fa)
subreads_bam=/home/shuaiw/borg/XRSBK_20221007_S64018_PL100268287-1_C01.subreads.bam
outdir=/home/shuaiw/borg/paper/ipdsummary/soil_1/

for ref in ${ref_set[@]}; do
    prefix=$outdir/$(basename $ref .fa)
    align_bam=$prefix.align.subreads.bam
    unsorted_bam=$prefix.align.unsorted.bam
    echo "    /usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref_dir/$ref $subreads_bam $unsorted_bam --preset SUBREAD -j 64"
    /usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref_dir/$ref $subreads_bam $unsorted_bam \
     --preset SUBREAD -j 64
    samtools sort -T $prefix -@ 64 -o "$align_bam" "$unsorted_bam"
    rm -f "$unsorted_bam"
    samtools index $align_bam
    /home/shuaiw//smrtlink/pbindex $align_bam
done

# $SLURM_CPUS_ON_NODE