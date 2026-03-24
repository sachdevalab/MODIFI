#!/bin/bash
#SBATCH --job-name=ipd3_test_400
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/ipdsummary/soil_subset/"
mkdir -p "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out3/"

/usr/bin/time -v -o "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out3//test_400.ipdSummary.time" \
    ~/smrtlink/ipdSummary "/home/shuaiw/borg/paper/ipdsummary/soil_subset//test_400.align.subreads.bam" \
    --reference "/home/shuaiw/borg/paper/ipdsummary/subset_ref//test_400.fa" \
    --numWorkers 64 \
    --gff "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out3//test_400.gff" \
    --csv "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out3//test_400.csv" \
    --methylFraction \
    --outfile "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out3//test_400.out"
