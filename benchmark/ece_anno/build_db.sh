#!/usr/bin/env bash
# Build the self-contained ECE-evidence HMM databases (run once).
#   - SCMG: gunzip + hmmpress bacterial (bac71) and archaeal (arc76) single-copy marker HMMs.
#   - Viral markers: extract a curated Pfam-A subset (terminase/capsid/portal) + hmmpress.
# Sources are copied into this repo dir so the module never depends on /home/mingy or /shared at run time.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# HMM DBs live in the data output dir (not the code dir). Override with $1.
DB_DIR="${1:-/home/shuaiw/borg/revision/ece_anno/db}"
SCMG_DIR="${DB_DIR}/scmg"
VIR_DIR="${DB_DIR}/viral_markers"
mkdir -p "${SCMG_DIR}" "${VIR_DIR}"

MING_DB="/home/mingy/MING/db"
PFAM_HMM="/shared/db/pfam/37/Pfam-A.hmm"

echo "=== [1/2] SCMG chromosomal single-copy marker HMMs ==="
for base in bac71 arc76; do
  src="${MING_DB}/${base}.hmm.gz"
  dst="${SCMG_DIR}/${base}.hmm"
  if [[ ! -f "$src" ]]; then echo "ERROR: missing $src" >&2; exit 1; fi
  echo "  gunzip ${src} -> ${dst}"
  zcat "$src" > "$dst"
  hmmpress -f "$dst" >/dev/null
  n=$(grep -c "^NAME " "$dst")
  echo "  ${base}: ${n} models, pressed."
done

echo "=== [2/2] Viral structural-marker Pfam subset ==="
VIR_HMM="${VIR_DIR}/viral_markers.hmm"
python "${SCRIPT_DIR}/build_viral_hmm.py" "${PFAM_HMM}" "${VIR_HMM}"
hmmpress -f "${VIR_HMM}" >/dev/null
echo "  viral_markers: $(grep -c '^NAME ' "${VIR_HMM}") models, pressed."

echo "=== done. DB layout: ==="
find "${DB_DIR}" -maxdepth 2 -type f | sort
