#!/bin/bash
#SBATCH --job-name=smk      # Job name
#SBATCH --partition=standard # Partition name

sbatch --partition standard --wrap "/usr/bin/time -v -o test.time snakemake "
