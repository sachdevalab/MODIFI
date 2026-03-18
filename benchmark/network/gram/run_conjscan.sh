#!/usr/bin/env bash
# Run ConjScan (macsyfinder with CONJScan models) on MGE protein FASTA.
# Usage: ./run_conjscan.sh [input.faa] [out_dir] [source_name] [threads] [models_dir]
# Default: input=.../gram_positive_linked_plasmids.faa_temp, out_dir=.../gram_positive, source_name=gram_positive_linked_plasmids, threads=48, models_dir=
# If models_dir is set, macsyfinder is called with --models-dir models_dir (path to CONJScan model package).
# Requires CONJScan models installed (e.g. macsydata install CONJScan) or pass models_dir. Output: best_solution_summary.tsv and all_systems.tsv.

set -euo pipefail

input_faa="${1:-/home/shuaiw/borg/paper/gram_positive/gram_positive_linked_plasmids.faa_temp}"
out_dir="${2:-/home/shuaiw/borg/paper/gram_positive}"
source_name="${3:-gram_positive_linked_plasmids}"
threads="${4:-48}"
models_dir="${5:-}"

params_out_dir="${out_dir}/${source_name}_anno/conjscan"
conjscan_summary="${params_out_dir}/best_solution_summary.tsv"
all_systems_tsv="${params_out_dir}/all_systems.tsv"

if [[ ! -s "${input_faa}" ]] || [[ "${source_name}" == "ncbi_chromosom_level_genomes" ]]; then
  echo "Skipping ConjScan: input empty or source_name is ncbi_chromosom_level_genomes."
  mkdir -p "${params_out_dir}"
  touch "${conjscan_summary}" "${all_systems_tsv}"
  exit 0
fi

MACSYFINDER="${MACSYFINDER:-/home/shuaiw/miniconda3/envs/conjscan/bin/macsyfinder}"
rm -rf "${params_out_dir}"
extra=()
[[ -n "${models_dir}" ]] && extra=(--models-dir "${models_dir}")
"${MACSYFINDER}" --db-type gembase -o "${params_out_dir}" --sequence-db "${input_faa}" -w "${threads}" "${extra[@]}" -m CONJScan
echo "ConjScan done. Summary: ${conjscan_summary}; all_systems: ${all_systems_tsv}"
