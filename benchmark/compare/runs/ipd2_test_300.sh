#!/bin/bash
#SBATCH --job-name=ipd2_test_300
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/ipdsummary/soil_1/"
mkdir -p "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2/"

/usr/bin/time -v -o "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2//test_300.ipdSummary.time" \
    ~/smrtlink/ipdSummary "/home/shuaiw/borg/paper/ipdsummary/soil_1//test_300.align.subreads.bam" \
    --reference "/home/shuaiw/borg/contigs//test_300.fa" \
    --numWorkers 64 \
    --gff "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2//test_300.gff" \
    --csv "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2//test_300.csv" \
    --methylFraction \
    --outfile "/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2//test_300.out"
