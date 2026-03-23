#!/usr/bin/env bash
# Regenerates ./runs/align2_*.sh and submits one Slurm job per reference.
# For submit-only: ./submit_align2_runs.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$SCRIPT_DIR/generate_align2_runs.sh"
"$SCRIPT_DIR/submit_align2_runs.sh"
