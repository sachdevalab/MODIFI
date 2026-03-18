#!/usr/bin/env bash
# Predict proteins from MGE/plasmid FASTA using pyrodigal-gv (meta mode), then sanitize headers.
# Usage: ./run_pyrodigal.sh [input.mge_fna] [threads]
# Default input: /home/shuaiw/borg/paper/gram_positive/gram_positive_linked_plasmids.fasta, threads=8

set -euo pipefail

input="${1:-/home/shuaiw/borg/paper/gram_positive/gram_positive_linked_plasmids.fasta}"
threads="${2:-8}"

base="${input%.*}"
output_faa="${base}.faa"
output_fna="${base}.fna"
output_gff="${base}.gff"
output_faa_temp="${base}.faa_temp"

if [[ ! -f "$input" ]]; then
  echo "Error: input not found: $input" >&2
  exit 1
fi

pyrodigal-gv -i "$input" -a "$output_faa" -d "$output_fna" -o "$output_gff" -p meta -j "$threads"
sed 's/#/_/g' "$output_faa" > "$output_faa_temp"
echo "Done. Proteins: $output_faa, sanitized: $output_faa_temp; genes: $output_fna; gff: $output_gff"
