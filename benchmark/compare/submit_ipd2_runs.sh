#!/usr/bin/env bash
# Submits all per-reference ipd2 jobs (run generate_ipd2_runs.sh first).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_DIR="$SCRIPT_DIR/runs"
ref_set=(test_100.fa test_200.fa test_300.fa test_500.fa)

missing=0
for ref in "${ref_set[@]}"; do
    base="$(basename "$ref" .fa)"
    f="$RUNS_DIR/ipd2_${base}.sh"
    if [[ ! -f "$f" ]]; then
        echo "Missing $f — run ./generate_ipd2_runs.sh first" >&2
        missing=1
    fi
done
if [[ "$missing" -ne 0 ]]; then
    exit 1
fi

for ref in "${ref_set[@]}"; do
    base="$(basename "$ref" .fa)"
    sbatch "$RUNS_DIR/ipd2_${base}.sh"
done
