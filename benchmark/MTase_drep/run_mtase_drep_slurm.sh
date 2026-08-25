#!/bin/bash
#SBATCH -p standard
#SBATCH -J mtase_drep
#SBATCH --qos unlimit-submit-20-run
#SBATCH -o /home/shuaiw/MODIFI/benchmark/MTase_drep/mtase_drep_%A_%a.out
#SBATCH --array=0-6

# Per-genome (per-contig) mmseqs dereplication of MTase genes + new MTase-vs-motif
# correlation, run as an array over amino-acid identity cutoffs (coverage fixed at
# 0.8) so we can compare and pick the best threshold.
# ~1180 tiny mmseqs easy-cluster jobs per cutoff, parallelized across contigs (each
# mmseqs job is --threads 1), so we use the whole node.
# NOTE: no subprocess/bash timeouts here; SLURM manages walltime.

set -euo pipefail

# identity cutoffs (fixed coverage 0.8); index by SLURM_ARRAY_TASK_ID
IDS=(0.30 0.50 0.70 0.80 0.90 0.95 1.00)
export MTASE_MIN_SEQ_ID="${IDS[$SLURM_ARRAY_TASK_ID]}"
export MTASE_MIN_COV="0.8"

# node-local scratch for the transient per-contig mmseqs churn (fast, not NFS)
export MTASE_DREP_WORK="/tmp/mtase_drep_work_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$MTASE_DREP_WORK"
trap 'rm -rf "$MTASE_DREP_WORK"' EXIT

echo "host=$(hostname) cpus=${SLURM_CPUS_ON_NODE} id=${MTASE_MIN_SEQ_ID} cov=${MTASE_MIN_COV} work=${MTASE_DREP_WORK}"
python -u /home/shuaiw/MODIFI/benchmark/MTase_drep/mtase_drep_corr.py
