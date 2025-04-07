#!/bin/bash

# HiFi Read Downsampling, Assembly, and Mapping Pipeline
# Author: YOU
# Usage: ./hifi_pipeline.sh input.bam output_prefix threads

# set -euo pipefail

# ----------------------
# Configurable Inputs
# ----------------------
# INPUT_BAM=$1            # Path to input BAM file
# PREFIX=$2               # Output prefix
# THREADS=${3:-16}        # Number of threads (default 16)

INPUT_BAM=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.hifi_reads.bam
PREFIX=/home/shuaiw/methylation/data/borg/bench/p5/p5
THREADS=10

# Derived filenames
DOWNSAMPLED_BAM="${PREFIX}.downsampled.bam"
FASTQ="${PREFIX}.fastq"
ASSEMBLY_DIR="${PREFIX}_assembly"
ASSEMBLY_FASTA="${ASSEMBLY_DIR}.p_ctg.fa"
MAPPED_BAM="${PREFIX}.mapped.bam"

# ----------------------
# Step 1: Downsample BAM to 5%
# ----------------------
# echo "[Step 1] Downsampling BAM to 5%..."
# samtools view -s 0.05 -b "$INPUT_BAM" -o "$DOWNSAMPLED_BAM"

# ----------------------
# Step 2: Convert BAM to FASTQ
# ----------------------
# echo "[Step 2] Converting downsampled BAM to FASTQ..."
# ## use samtools fastq
# samtools fastq "$DOWNSAMPLED_BAM" > "$FASTQ"
# gzip "$FASTQ" # Compress the FASTQ file
# # Result: $FASTQ or $PREFIX.fastq

# ----------------------
# Step 3: Assemble with hifiasm
# ----------------------
echo "[Step 3] Assembling reads with hifiasm..."
mkdir -p "$ASSEMBLY_DIR"
hifiasm_meta -o "$ASSEMBLY_DIR/$PREFIX" -t "$THREADS" "$FASTQ.gz" > "$ASSEMBLY_DIR/asm.log"

# ----------------------
# Step 4: Map Reads Back to Assembly (using pbmm2)
# ----------------------
echo "[Step 4] Mapping reads back to contigs with pbmm2..."
~/smrtlink/pbmm2 align "$ASSEMBLY_FASTA" "$DOWNSAMPLED_BAM" "$MAPPED_BAM" \
    --preset CCS --sort --threads "$THREADS"
samtools index "$MAPPED_BAM"
/home/shuaiw//smrtlink/pbindex "$MAPPED_BAM"

echo "[Done] Results:"
echo "  - Downsampled BAM: $DOWNSAMPLED_BAM"
echo "  - FASTQ: $FASTQ"
echo "  - Assembly FASTA: $ASSEMBLY_FASTA"
echo "  - Mapped BAM: $MAPPED_BAM"

