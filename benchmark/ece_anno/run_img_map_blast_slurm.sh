#!/bin/bash
#SBATCH -p standard
#SBATCH -J img_map_blast
#SBATCH -o /home/shuaiw/borg/revision/ece_anno/expanded/img_map/blast_%j.out
#SBATCH --qos unlimit-submit-20-run

set -euo pipefail
D=/home/shuaiw/borg/revision/ece_anno/expanded/img_map
IMGVR=/shared/db/imgvr/v4_IMG_VR_2022-09-20_6/nucleotide/IMGVR_all_nucleotides.fna
IMGPR=/shared/db/imgpr/2023-08-08_1/IMGPR_nucl.fna
BL=/shared/software/bin/blastn
T=$SLURM_CPUS_ON_NODE

echo "[`date`] viruses vs IMG/VR ($T threads)"
$BL -task megablast -query $D/linked_viruses.fna -db $IMGVR \
  -outfmt '6 std qlen slen' -evalue 1e-10 -max_target_seqs 25 -num_threads $T \
  > $D/virus_imgvr_blast.tsv
echo "[`date`] plasmids vs IMG/PR ($T threads)"
$BL -task megablast -query $D/linked_plasmids.fna -db $IMGPR \
  -outfmt '6 std qlen slen' -evalue 1e-10 -max_target_seqs 25 -num_threads $T \
  > $D/plasmid_imgpr_blast.tsv
echo "[`date`] done. virus hits=$(wc -l < $D/virus_imgvr_blast.tsv) plasmid hits=$(wc -l < $D/plasmid_imgpr_blast.tsv)"

# ANI (cluster_MGE.py stringency) -- anicalc computes weighted ANI (pid) + qcov + tcov per (query,ref) pair
PY=/home/shuaiw/miniconda3/envs/methy3/bin/python
ANICALC=/groups/sachdeva/projects/sag/SAGLink/workflow/anicalc.py
$PY $ANICALC -i $D/virus_imgvr_blast.tsv   -o $D/virus_imgvr_ani.tsv
$PY $ANICALC -i $D/plasmid_imgpr_blast.tsv -o $D/plasmid_imgpr_ani.tsv
echo "[`date`] anicalc done. virus_ani=$(wc -l < $D/virus_imgvr_ani.tsv) plasmid_ani=$(wc -l < $D/plasmid_imgpr_ani.tsv)"
