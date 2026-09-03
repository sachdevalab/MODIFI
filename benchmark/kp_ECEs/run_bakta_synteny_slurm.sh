#!/bin/bash
#SBATCH -p standard
#SBATCH -J bakta_syn
#SBATCH -o /home/shuaiw/borg/revision/kp_eces/synteny/bakta_syn_%j.out
#SBATCH --qos unlimit-submit-20-run
# Annotate the 18 infant_15_35_C cluster members with bakta (light DB v6, compatible with bakta 1.12.1).
# Full-node allocation: run several members in parallel, each with a share of the node's cores.
set -uo pipefail
export PATH=/shared/software/bin:$PATH   # compute node lacks some tool paths (pilercr etc.)
OUT=/home/shuaiw/borg/revision/kp_eces/synteny
DB=/home/allisong/bin/bakta_db/db-light
mkdir -p "$OUT/bakta" "$OUT/gbff"

NCPU=${SLURM_CPUS_ON_NODE:-32}
JOBS=6                         # 6 members at a time
THREADS=$(( NCPU / JOBS )); [ "$THREADS" -lt 1 ] && THREADS=1
echo "cores=$NCPU jobs=$JOBS threads/job=$THREADS"

run_one() {
  m=$1; DB=$2; OUT=$3; THREADS=$4
  lt=$(echo "$m" | tr -d '_')
  conda run -n bakta_env bakta --db "$DB" --keep-contig-headers --skip-crispr \
    --output "$OUT/bakta/$m" --prefix "$m" --locus-tag "$lt" \
    --threads "$THREADS" --force "$OUT/contigs/$m.fna" \
    > "$OUT/bakta/$m.log" 2>&1 && echo "OK $m" || echo "FAIL $m"
}
export -f run_one

ls "$OUT"/contigs/*.fna | sed 's|.*/||; s|\.fna$||' | \
  xargs -I{} -P "$JOBS" bash -c 'run_one "$@"' _ {} "$DB" "$OUT" "$THREADS"

# collect gbff
for f in "$OUT"/bakta/*/*.gbff; do cp "$f" "$OUT/gbff/"; done
echo "gbff collected: $(ls $OUT/gbff/*.gbff 2>/dev/null | wc -l) / 18"
