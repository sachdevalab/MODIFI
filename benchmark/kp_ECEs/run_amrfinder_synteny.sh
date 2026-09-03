#!/bin/bash
#SBATCH -p standard
#SBATCH -J amr_syn
#SBATCH -o /home/shuaiw/borg/revision/kp_eces/synteny/amr_syn_%j.out
#SBATCH --qos unlimit-submit-20-run
# AMRFinderPlus (nucleotide) on the 18 members combined -> contig coordinates of AMR + STRESS genes.
set -uo pipefail
OUT=/home/shuaiw/borg/revision/kp_eces/synteny
NCPU=${SLURM_CPUS_ON_NODE:-16}
conda run -n kp_eces amrfinder -n "$OUT/members.fna" --plus --threads "$NCPU" \
  -o "$OUT/synteny_amrfinder.tsv" > "$OUT/amrfinder.log" 2>&1
echo "exit=$?"
echo "rows: $(($(wc -l < $OUT/synteny_amrfinder.tsv) - 1))"
cut -f10 "$OUT/synteny_amrfinder.tsv" | tail -n +2 | sort | uniq -c
