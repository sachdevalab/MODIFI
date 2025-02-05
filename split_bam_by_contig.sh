#!/bin/bash

# Input BAM file
input_bam="/home/shuaiw/borg/break_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam"
output_dir="/home/shuaiw/methylation/data/borg/b_contigs/bams"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Get the list of contigs
contigs=$(samtools idxstats "$input_bam" | cut -f 1 | grep -v '*')

# Loop through each contig and create a separate BAM file
for contig in $contigs; do
    output_bam="$output_dir/${contig}.bam"
    echo "Processing contig: $contig"
    samtools view -b "$input_bam" "$contig" > "$output_bam"
    samtools index "$output_bam"
done

echo "Splitting completed."