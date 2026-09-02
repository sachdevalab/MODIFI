#!/bin/bash
#SBATCH --job-name=ft_scale
#SBATCH --partition=memory
#SBATCH --qos=unlimit-submit-20-run
#SBATCH --array=0-2
#SBATCH -o /home/shuaiw/borg/paper/ipdsummary/soil_1/fibertools.out/slurm_%A_%a.out
#SBATCH -e /home/shuaiw/borg/paper/ipdsummary/soil_1/fibertools.out/slurm_%A_%a.err
# Run the fibertools scaling points test_200/300/500 (one per array task) on the
# memory partition (64-core nodes). THREADS is fixed to 64 inside batch_fibertools.sh
# to match the test_100 terminal run so wall-clock times stay comparable.
REFS=(test_200 test_300 test_500)
base=${REFS[$SLURM_ARRAY_TASK_ID]}
echo "[ft_scale] $base on $(hostname) cpus=$SLURM_CPUS_ON_NODE $(date)"
bash /home/shuaiw/MODIFI/benchmark/compare/batch_fibertools.sh "$base"
echo "[ft_scale] $base finished $(date)"
