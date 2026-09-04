#!/bin/bash
# antiSMASH 7.1.0 BGC detection on the 14 strict-set Kp ECE representatives (terminal, not SLURM).
# Per-rep run for clean attribution; a few in parallel with modest cpus each.
set -uo pipefail
export PATH=/shared/software/bin:$PATH
GP=/home/shuaiw/borg/revision/kp_eces/gene_profile
OUT=/home/shuaiw/borg/revision/kp_eces/antismash
mkdir -p "$OUT/contigs"

# ensure per-rep FASTAs exist
for r in $(awk -F'\t' 'NR>1{print $2}' "$GP/kp14_representatives.tsv"); do
  [ -s "$OUT/contigs/$r.fna" ] || conda run -n tldr seqkit grep -p "$r" "$GP/kp14_reps.fna" > "$OUT/contigs/$r.fna"
done

run_one() {
  r=$1; OUT=$2
  [ -s "$OUT/$r/$r.json" ] && { echo "SKIP $r (done)"; return; }
  antismash --taxon bacteria --genefinding-tool prodigal --cb-knownclusters --cb-general \
    --output-dir "$OUT/$r" --cpus 8 "$OUT/contigs/$r.fna" > "$OUT/$r.log" 2>&1 \
    && echo "OK $r" || echo "FAIL $r"
}
export -f run_one

awk -F'\t' 'NR>1{print $2}' "$GP/kp14_representatives.tsv" | \
  xargs -I{} -P 4 bash -c 'run_one "$@"' _ {} "$OUT"
echo "antiSMASH done: $(ls -d $OUT/*/ 2>/dev/null | wc -l) output dirs"
