#!/bin/bash
#SBATCH --job-name=align2_test_100
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/ipdsummary/soil_subset"

/usr/bin/time -v -o "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100.pbmm2.align.time" \
    ~/smrtlink/pbmm2 align "/home/shuaiw/borg/paper/ipdsummary/subset_ref/test_100.fa" "/home/shuaiw/borg/XRSBK_20221007_S64018_PL100268287-1_C01.subreads.bam" "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100.align.unsorted.bam" \
    --preset SUBREAD -j 64
samtools sort -T "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100" -@ 64 -o "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100.align.subreads.bam" "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100.align.unsorted.bam"
rm -f "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100.align.unsorted.bam"
samtools index "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100.align.subreads.bam"
~/smrtlink/pbindex "/home/shuaiw/borg/paper/ipdsummary/soil_subset/test_100.align.subreads.bam"
