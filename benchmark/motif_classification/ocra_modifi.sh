#!/bin/bash
#SBATCH -p standard
#SBATCH -J ocra_modifi
#SBATCH --qos unlimit-submit-20-run
#SBATCH -o /home/shuaiw/borg/revision/ocra_5mC/ocra_modifi_%x_%j.slurm.out
# Usage: sbatch ocra_modifi.sh <input_subreads.bam> <outdir_name>
set -euo pipefail
source /home/shuaiw/miniconda3/etc/profile.d/conda.sh
conda activate MODIFI_subreads
export PATH=/home/shuaiw/miniconda3/envs/MODIFI_subreads/bin:$PATH
R=/home/shuaiw/borg/revision/ocra_5mC
INBAM="$1"; OUTNAME="$2"
KMEAN=/home/shuaiw/borg/paper/gg_run3/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat
KNUM=/home/shuaiw/borg/paper/gg_run3/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat

# subreads mode; soil_1 external k-mer control (isolate data too simple for a self control);
# permissive motif thresholds so a weak 5mC motif is not filtered out before inspection.
python /home/shuaiw/MODIFI/main.py \
  -b "$R/$INBAM" \
  -r "$R/ocra_K2.fa" \
  -o "$R/$OUTNAME/" \
  --read_type subreads \
  --min_frac 0.2 --min_sites 50 \
  --kmer_mean_db "$KMEAN" \
  --kmer_num_db  "$KNUM" \
  --threads $SLURM_CPUS_ON_NODE
echo "MODIFI DONE: $OUTNAME"
