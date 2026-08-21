#!/bin/bash
#SBATCH -p standard
#SBATCH -J kp_deepdive
#SBATCH -o /home/shuaiw/borg/revision/kp_eces/slurm_out/deepdive_%j.out
# Step 3 - deep dive on the infant_15_35_C plasmid cluster (reviewer part 3).
# Annotates the 416 kb representative (Citrobacter host) and the Klebsiella member infant_8_26_C
# (112 kb) with bakta (genes/products/AMR/replicons) + mob_typer (mobility, replicon/relaxase MOB
# type, closest known plasmid + mash distance + predicted host range), and builds the cross-genus
# member table (each member's host species/genus/length + ANI to the representative).
# Whole node, no subprocess timeouts.
set -uo pipefail

source /shared/software/miniconda3/latest/etc/profile.d/conda.sh
conda activate kp_eces
ENV=/home/shuaiw/miniconda3/envs/kp_eces
BENV=/home/shuaiw/miniconda3/envs/bakta_env
DBROOT=/home/shuaiw/borg/revision/kp_eces/dbs
BAKTA_DB=/home/shuaiw/borg/revision/kp_eces/dbs/bakta/db-light
SEQKIT=/home/shuaiw/miniconda3/envs/tldr/bin/seqkit
ALLMGE=/home/shuaiw/borg/paper/network/all_mge.fa
ANIFILE=/home/shuaiw/borg/paper/MGE/cluster/megablast.ani.95ani.tsv
OUT=/home/shuaiw/borg/revision/kp_eces/deepdive
FIG=/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs
mkdir -p "$OUT"
THREADS=${SLURM_CPUS_ON_NODE:-16}
export PATH="$ENV/bin:$PATH"

for seq in infant_15_35_C infant_8_26_C; do
  echo "[deepdive] $(date) === $seq ==="
  $SEQKIT grep -p "$seq" "$ALLMGE" > "$OUT/$seq.fna" 2>/dev/null
  echo "[deepdive] length: $($SEQKIT stats -T "$OUT/$seq.fna" 2>/dev/null | tail -1 | cut -f5)"

  # mob_typer
  "$ENV/bin/mob_typer" --infile "$OUT/$seq.fna" --out_file "$OUT/$seq.mobtyper.txt" \
      --database_directory "$DBROOT/mob_suite" 2> "$OUT/$seq.mobtyper.log" || \
      echo "[deepdive] mob_typer failed for $seq"

  # bakta full annotation (light DB), only if the DB finished downloading
  if [ -x "$BENV/bin/bakta" ] && [ -d "$BAKTA_DB" ]; then
    "$BENV/bin/bakta" --db "$BAKTA_DB" --skip-plot --force \
        --output "$OUT/${seq}_bakta" --prefix "$seq" --threads "$THREADS" \
        "$OUT/$seq.fna" 2> "$OUT/$seq.bakta.log" || echo "[deepdive] bakta failed for $seq"
  else
    echo "[deepdive] bakta DB not ready ($BAKTA_DB); skipping full gene annotation for $seq"
  fi
done

# cross-genus member table: ANI of each infant_15_35_C member to the representative
echo "[deepdive] $(date) building infant_15_35_C member ANI table"
awk -F'\t' 'BEGIN{OFS="\t"} $1=="infant_15_35_C" || $2=="infant_15_35_C"' "$ANIFILE" \
    > "$OUT/infant_15_35_C.ani_rows.tsv" 2>/dev/null || true
echo "[deepdive] ANI rows for representative: $(wc -l < "$OUT/infant_15_35_C.ani_rows.tsv" 2>/dev/null || echo 0)"

echo "[deepdive] $(date) DONE"
