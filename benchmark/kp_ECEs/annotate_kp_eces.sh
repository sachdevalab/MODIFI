#!/bin/bash
#SBATCH -p standard
#SBATCH -J kp_annotate
#SBATCH -o /home/shuaiw/borg/revision/kp_eces/slurm_out/annotate_%j.out
# Step 2 - known-plasmid identity, mobility, host range, and AMR for the Klebsiella-associated ECEs.
# mob_typer (--multi: mobility + replicon/relaxase MOB type + closest reference plasmid + mash
# distance + predicted host range) on all cluster-member ECEs, and AMRFinderPlus (AMR + virulence +
# stress) on the same set. Whole node ($SLURM_CPUS_ON_NODE). NO subprocess timeouts.
set -uo pipefail

source /shared/software/miniconda3/latest/etc/profile.d/conda.sh
conda activate kp_eces
DBROOT=/home/shuaiw/borg/revision/kp_eces/dbs
SEQS=/home/shuaiw/borg/revision/kp_eces/seqs
OUT=/home/shuaiw/borg/revision/kp_eces/annotate
mkdir -p "$OUT"
THREADS=${SLURM_CPUS_ON_NODE:-16}

MODE="${1:-all}"   # all | mob | amr

if [ "$MODE" = "all" ] || [ "$MODE" = "mob" ]; then
  echo "[annotate] $(date) mob_typer --multi on cluster members ($THREADS threads)"
  mob_typer --multi --infile "$SEQS/kp_cluster_members.fna" \
      --out_file "$OUT/mob_typer.tsv" --num_threads "$THREADS" \
      --database_directory "$DBROOT/mob_suite" 2> "$OUT/mob_typer.log"
  echo "[annotate] mob_typer rows: $(($(wc -l < "$OUT/mob_typer.tsv" 2>/dev/null || echo 1) - 1))"
fi

if [ "$MODE" = "all" ] || [ "$MODE" = "amr" ]; then
  echo "[annotate] $(date) AMRFinderPlus (nucleotide mode, --plus) on cluster members"
  # -u downloads to the env default DB dir; do NOT pass -d (rejected with -u). Use default here too.
  amrfinder -n "$SEQS/kp_cluster_members.fna" --plus --threads "$THREADS" \
      --name kp_eces -o "$OUT/amrfinder.tsv" 2> "$OUT/amrfinder.log"
  echo "[annotate] AMR hits: $(($(wc -l < "$OUT/amrfinder.tsv" 2>/dev/null || echo 1) - 1))"
fi

echo "[annotate] $(date) DONE ($MODE)"
