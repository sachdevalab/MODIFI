#!/bin/bash
#SBATCH --job-name=ipd_summary_28   # Job name
#SBATCH --partition=standard        # Partition name

set -euo pipefail

ref_dir=/home/shuaiw/borg/contigs/
ref_set=(test_100.fa test_200.fa test_300.fa test_500.fa)
outdir=/home/shuaiw/borg/paper/ipdsummary/soil_1/
ipd_output=/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2/

mkdir -p "$outdir"
mkdir -p "$ipd_output"

for ref in "${ref_set[@]}"; do
    prefix="$outdir/$(basename "$ref" .fa)"
    ipd_prefix="$ipd_output/$(basename "$ref" .fa)"
    align_bam="$prefix.align.subreads.bam"
    ref_path="$ref_dir/$ref"

    echo "Running ipdSummary for $ref_path"
    echo "    /usr/bin/time -v -o $ipd_prefix.ipdSummary.time ~/smrtlink/ipdSummary $align_bam --reference $ref_path --debug --numWorkers 1 --gff $ipd_prefix.gff --csv $ipd_prefix.csv --methylFraction --outfile $ipd_prefix.out"

    /usr/bin/time -v -o "$ipd_prefix.ipdSummary.time" \
        ~/smrtlink/ipdSummary "$align_bam" \
        --reference "$ref_path" \
        --numWorkers 64 \
        --gff "$ipd_prefix.gff" \
        --csv "$ipd_prefix.csv" \
        --methylFraction \
        --outfile "$ipd_prefix.out"

    # break
done
