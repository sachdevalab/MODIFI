#!/bin/bash
#SBATCH -p standard
#SBATCH -J jf8_wga
#SBATCH --qos unlimit-submit-20-run
#SBATCH -o /home/shuaiw/borg/revision/motif_class/jf8_wga.slurm.out

# JF8 mock rerun, subreads mode, using the 3-strain RS II WGA k-mer control DB
# (chemistry-matched, modification-free). New output dir; existing jf8_self/jf8_soil1 untouched.
source /home/shuaiw/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate MODIFI_subreads 2>/dev/null
export PATH=/home/shuaiw/miniconda3/envs/MODIFI_subreads/bin:$PATH

python /home/shuaiw/MODIFI/main.py \
  --aligned_bam /home/shuaiw/methylation/data/published_data/fanggang/align/Mock_JF8.align.bam \
  -r /home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa \
  -o /home/shuaiw/borg/revision/motif_class/jf8_wga/ \
  --read_type subreads \
  --threads $SLURM_CPUS_ON_NODE \
  --kmer_mean_db /home/shuaiw/borg/revision/ocra_5mC/wga_control/control/control_db.up7.down3.mean.dat \
  --kmer_num_db  /home/shuaiw/borg/revision/ocra_5mC/wga_control/control/control_db.up7.down3.num.dat
