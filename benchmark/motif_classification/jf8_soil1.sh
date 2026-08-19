#!/bin/bash
#SBATCH -p standard
#SBATCH -J jf8_soil1
#SBATCH --qos unlimit-submit-20-run
#SBATCH -o /home/shuaiw/borg/revision/motif_class/jf8_soil1.slurm.out

# Use the MODIFI_subreads conda env (has adjustText, pbmm2, pbindex).
# pyenv shims shadow conda's python on PATH, so prepend the env bin explicitly.
source /home/shuaiw/miniconda3/etc/profile.d/conda.sh
conda activate MODIFI_subreads
export PATH=/home/shuaiw/miniconda3/envs/MODIFI_subreads/bin:$PATH

# JF8 mock rerun, subreads mode, SOIL_1 k-mer control database (high-complexity background;
# JF8's own diversity is too low to build a good self control model)
python /home/shuaiw/MODIFI/main.py \
  --aligned_bam /home/shuaiw/methylation/data/published_data/fanggang/align/Mock_JF8.align.bam \
  -r /home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa \
  -o /home/shuaiw/borg/revision/motif_class/jf8_soil1/ \
  --read_type subreads \
  --threads $SLURM_CPUS_ON_NODE \
  --kmer_mean_db /home/shuaiw/borg/paper/gg_run3/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
  --kmer_num_db  /home/shuaiw/borg/paper/gg_run3/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat
