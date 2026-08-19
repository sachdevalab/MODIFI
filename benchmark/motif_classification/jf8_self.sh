#!/bin/bash
#SBATCH -p standard
#SBATCH -J jf8_self
#SBATCH --qos unlimit-submit-20-run
#SBATCH -o /home/shuaiw/borg/revision/motif_class/jf8_self.slurm.out

# Use the MODIFI_subreads conda env (has adjustText, pbmm2, pbindex).
# pyenv shims shadow conda's python on PATH, so prepend the env bin explicitly.
source /home/shuaiw/miniconda3/etc/profile.d/conda.sh
conda activate MODIFI_subreads
export PATH=/home/shuaiw/miniconda3/envs/MODIFI_subreads/bin:$PATH

# JF8 mock rerun, subreads mode, SELF k-mer control database (built from JF8 itself)
python /home/shuaiw/MODIFI/main.py \
  --aligned_bam /home/shuaiw/methylation/data/published_data/fanggang/align/Mock_JF8.align.bam \
  -r /home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa \
  -o /home/shuaiw/borg/revision/motif_class/jf8_self/ \
  --read_type subreads \
  --threads $SLURM_CPUS_ON_NODE
