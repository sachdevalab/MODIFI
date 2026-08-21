#!/bin/bash
#SBATCH -p standard
#SBATCH -J phasevar_orient
#SBATCH -o /home/shuaiw/borg/revision/phase_variation/orientation/orient_%j.out

# Read-level orientation of the invertible TRD across time points.
# Usage: sbatch run_orientation_slurm.sh [cluster]   (omit cluster to run all 6)
# SLURM manages walltime via -t / full node; NO subprocess timeouts here.

set -euo pipefail
source /shared/software/miniconda3/latest/etc/profile.d/conda.sh
conda activate modifi

export PATH=/shared/software/bin:$PATH   # minimap2, samtools, seqkit
mkdir -p /home/shuaiw/borg/revision/phase_variation/orientation

cd /home/shuaiw/MODIFI/benchmark/phase_variation
echo "host=$(hostname) cpus=${SLURM_CPUS_ON_NODE:-NA} cluster_arg=${1:-ALL}"
python orientation_driver.py ${1:-}
echo "DONE"
