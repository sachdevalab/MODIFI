#!/bin/bash
#SBATCH -p standard
#SBATCH -J hostmotif_density
#SBATCH -o /home/shuaiw/MODIFI/tmp/rev_figs/borg/hostmotif_density_%j.out
#SBATCH --qos unlimit-submit-20-run

cd /home/shuaiw/MODIFI/benchmark/borg

/home/shuaiw/miniconda3/envs/modifi/bin/python motif_density_hostmotifs_per_sample.py \
    --threads "$SLURM_CPUS_ON_NODE"
