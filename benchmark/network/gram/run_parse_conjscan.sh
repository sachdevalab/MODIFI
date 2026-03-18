#!/usr/bin/env bash
# Parse ConjScan results (best_solution_summary.tsv + all_systems.tsv) into conjscan_parsed and conjscan_hits.
# Usage: ./run_parse_conjscan.sh [out_dir] [source_name]
# Default: out_dir=/home/shuaiw/borg/paper/gram_positive, source_name=gram_positive_linked_plasmids
# Expects ConjScan to have been run so that:
#   {out_dir}/{source_name}_anno/conjscan/best_solution_summary.tsv
#   {out_dir}/{source_name}_anno/conjscan/all_systems.tsv
# exist. Outputs:
#   {out_dir}/{source_name}_anno/conjscan/{source_name}_MGE_contigs.conjscan_parsed.tsv
#   {out_dir}/{source_name}_anno/conjscan/{source_name}_MGE_contigs.conjscan_hits.tsv

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out_dir="${1:-/home/shuaiw/borg/paper/gram_positive}"
source_name="${2:-gram_positive_linked_plasmids}"

conjscan_dir="${out_dir}/${source_name}_anno/conjscan"
conjscan_summary="${conjscan_dir}/best_solution_summary.tsv"
all_systems_tsv="${conjscan_dir}/all_systems.tsv"
conjscan_parsed="${conjscan_dir}/${source_name}_MGE_contigs.conjscan_parsed.tsv"
conjscan_hits="${conjscan_dir}/${source_name}_MGE_contigs.conjscan_hits.tsv"

if [[ ! -s "${conjscan_summary}" ]]; then
  echo "ConjScan summary missing or empty: ${conjscan_summary}; touching outputs and exiting."
  mkdir -p "${conjscan_dir}"
  touch "${conjscan_parsed}" "${conjscan_hits}"
  exit 0
fi

python /home/mingy/MING/helper_scripts/parse_conjscan.py \
  "${conjscan_summary}" \
  "${conjscan_parsed}" \
  "${all_systems_tsv}" \
  "${conjscan_hits}"
