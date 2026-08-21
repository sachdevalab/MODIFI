#!/bin/bash
# mob_typer (--multi) + AMRFinderPlus on all linked ECEs. Run in terminal (only 379).
set -uo pipefail
source /shared/software/miniconda3/latest/etc/profile.d/conda.sh
conda activate kp_eces
DBROOT=/home/shuaiw/borg/revision/kp_eces/dbs
SEQS=/home/shuaiw/borg/revision/linked_eces/seqs
OUT=/home/shuaiw/borg/revision/linked_eces/annotate
mkdir -p "$OUT"
THREADS=$(nproc)

echo "[linked-annotate] $(date) mob_typer --multi ($THREADS threads)"
mob_typer --multi --infile "$SEQS/linked_eces.fna" --out_file "$OUT/mob_typer.tsv" \
    --num_threads "$THREADS" --database_directory "$DBROOT/mob_suite" 2> "$OUT/mob_typer.log"
echo "[linked-annotate] mob_typer rows: $(($(wc -l < "$OUT/mob_typer.tsv") - 1))"

echo "[linked-annotate] $(date) AMRFinderPlus (nucleotide, --plus)"
amrfinder -n "$SEQS/linked_eces.fna" --plus --threads "$THREADS" \
    --name linked -o "$OUT/amrfinder.tsv" 2> "$OUT/amrfinder.log"
echo "[linked-annotate] AMR hits: $(($(wc -l < "$OUT/amrfinder.tsv") - 1))"
echo "[linked-annotate] $(date) DONE"
