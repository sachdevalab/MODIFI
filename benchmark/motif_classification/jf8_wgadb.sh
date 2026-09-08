#!/bin/bash
#SBATCH -p standard
#SBATCH -J jf8_wgadb
#SBATCH --qos unlimit-submit-20-run
#SBATCH --exclude=node-48-256g-12,node-48-256g-13,node-48-256g-14,node-48-256g-17,node-48-256g-18,node-48-256g-19,node-48-256g-2,node-48-256g-20,node-48-256g-3,node-48-256g-8,node-48-256g-9,node-48-384g-1,node-48-384g-2,node-48-384g-3,node-48-384g-4
#SBATCH -o /home/shuaiw/borg/revision/motif_class/jf8_wgadb.slurm.out
source /home/shuaiw/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate MODIFI_subreads 2>/dev/null
export PATH=/home/shuaiw/miniconda3/envs/MODIFI_subreads/bin:$PATH
python /home/shuaiw/MODIFI/main.py \
  --aligned_bam /home/shuaiw/methylation/data/published_data/fanggang/align/Mock_JF8.align.bam \
  -r /home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa \
  -o /home/shuaiw/borg/revision/motif_class/jf8_wgadb/ \
  --read_type subreads --no-clean \
  --kmer_mean_db /home/shuaiw/MODIFI/control_db/control_db.RSII.up7.down3.mean.dat --kmer_num_db /home/shuaiw/MODIFI/control_db/control_db.RSII.up7.down3.num.dat \
  --threads $SLURM_CPUS_ON_NODE
