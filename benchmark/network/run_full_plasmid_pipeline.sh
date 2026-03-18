#!/usr/bin/env bash
# Run full plasmid pipeline for a domain: classify (FASTA) -> pyrodigal -> conjscan -> parse_conjscan.
# Usage: ./run_full_plasmid_pipeline.sh <domain_suffix> [threads]
#   domain_suffix: gram_positive | gram_negative | archaea
#   threads: optional, default 48 (used for pyrodigal and conjscan)
# Example: ./run_full_plasmid_pipeline.sh gram_negative 48

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
domain_suffix="${1:?Usage: $0 <gram_positive|gram_negative|archaea> [threads]}"
threads="${2:-48}"

BORG_PAPER="/home/shuaiw/borg/paper"
case "${domain_suffix}" in
  gram_positive)  domain="Gram-Positive" ;;
  gram_negative)  domain="Gram-Negative" ;;
  archaea)        domain="Archaea" ;;
  *) echo "Error: domain_suffix must be gram_positive, gram_negative, or archaea" >&2; exit 1 ;;
esac

out_dir="${BORG_PAPER}/${domain_suffix}"
source_name="${domain_suffix}_linked_plasmids"
fasta="${out_dir}/${source_name}.fasta"
faa_temp="${out_dir}/${source_name}.faa_temp"

echo "=== Plasmid pipeline for ${domain} (${domain_suffix}) ==="
# 1. Collect FASTA
python "${SCRIPT_DIR}/classify_plasmids.py" --domain "${domain}"
if [[ ! -s "${fasta}" ]]; then
  echo "No FASTA produced (empty or missing). Skipping pyrodigal/conjscan/parse."
  exit 0
fi
# 2. Pyrodigal
"${SCRIPT_DIR}/gram/run_pyrodigal.sh" "${fasta}" "${threads}"
# 3. ConjScan
"${SCRIPT_DIR}/gram/run_conjscan.sh" "${faa_temp}" "${out_dir}" "${source_name}" "${threads}"
# 4. Parse ConjScan
"${SCRIPT_DIR}/gram/run_parse_conjscan.sh" "${out_dir}" "${source_name}"
echo "=== Done: ${domain_suffix} ==="
