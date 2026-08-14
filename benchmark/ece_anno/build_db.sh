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
PLAS_DIR="${DB_DIR}/plasmid_markers"
RRNA_DIR="${DB_DIR}/rrna"
mkdir -p "${SCMG_DIR}" "${VIR_DIR}" "${PLAS_DIR}" "${RRNA_DIR}"

MING_DB="/home/mingy/MING/db"
PFAM_HMM="/shared/db/pfam/37/Pfam-A.hmm"
RFAM_CM="/shared/db/rfam/latest/cm/Rfam.cm"
# Bacterial + archaeal SSU (16S) and LSU (23S) rRNA Rfam families
RRNA_RFAM="RF00177 RF02541 RF01959 RF02540"

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

echo "=== [2/4] Viral structural-marker Pfam subset ==="
VIR_HMM="${VIR_DIR}/viral_markers.hmm"
python "${SCRIPT_DIR}/build_viral_hmm.py" "${PFAM_HMM}" "${VIR_HMM}" viral
hmmpress -f "${VIR_HMM}" >/dev/null
echo "  viral_markers: $(grep -c '^NAME ' "${VIR_HMM}") models, pressed."

echo "=== [2b] VOGdb -> structural-class map (terminase/capsid/portal) ==="
VOG_ANN="/shared/db/vogdb/r98/metadata/vog.annotations.tsv"
python - "$VOG_ANN" "${VIR_DIR}/vog_structural_map.tsv" <<'PYEOF'
import pandas as pd, re, sys
a = pd.read_csv(sys.argv[1], sep="\t"); a.columns = [c.lstrip("#") for c in a.columns]
pats = {"terminase_large": r"terminase.*large|large.*terminase|\bTerL\b",
        "terminase_small": r"terminase.*small|small.*terminase|\bTerS\b",
        "major_capsid": r"capsid|major head|prohead", "portal": r"portal"}
rows = []
for _, r in a.iterrows():
    d = str(r["ConsensusFunctionalDescription"])
    for c, p in pats.items():
        if re.search(p, d, re.I):
            rows.append((r["GroupName"], c)); break
pd.DataFrame(rows, columns=["vog", "class"]).to_csv(sys.argv[2], sep="\t", index=False)
print(f"  vog_structural_map: {len(rows)} VOG->class")
PYEOF

echo "=== [3/4] Plasmid hallmark Pfam subset (Rep + ParA/ParB) ==="
PLAS_HMM="${PLAS_DIR}/plasmid_markers.hmm"
python "${SCRIPT_DIR}/build_viral_hmm.py" "${PFAM_HMM}" "${PLAS_HMM}" plasmid
hmmpress -f "${PLAS_HMM}" >/dev/null
echo "  plasmid_markers: $(grep -c '^NAME ' "${PLAS_HMM}") models, pressed."

echo "=== [4/4] rRNA Rfam covariance models (16S/23S, bact+arch) ==="
RRNA_CM="${RRNA_DIR}/rrna.cm"
if [[ ! -s "${RRNA_CM}" ]]; then
  : > "${RRNA_CM}"
  for rf in ${RRNA_RFAM}; do
    cmfetch "${RFAM_CM}" "${rf}" >> "${RRNA_CM}"
  done
  cmpress -F "${RRNA_CM}" >/dev/null
fi
echo "  rrna: $(grep -c '^NAME ' "${RRNA_CM}") models, pressed."

echo "=== done. DB layout: ==="
find "${DB_DIR}" -maxdepth 2 -type f | sort
