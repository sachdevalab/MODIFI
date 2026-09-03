#!/bin/bash
# Two-pass LoVis4u synteny plot of the 18 infant_15_35_C cluster members, resistance genes marked.
# Pass 1 (no mmseqs) -> exact feature_id + coordinates LoVis4u assigns; build_faf.py maps AMR/STRESS
# hits to those ids; Pass 2 draws the synteny with the -faf resistance marks.
set -euo pipefail
OUT=/home/shuaiw/borg/revision/kp_eces/synteny
FIG=/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs
CFG=/home/shuaiw/MODIFI/benchmark/kp_ECEs/lovis4u_synteny.cfg
GB=$OUT/gbff
L4="conda run -n lovis4u lovis4u"
mkdir -p "$FIG"

echo "[clean] drop redundant bakta 'gene' features (duplicate locus_tags) -> gbff_clean"
CLEAN=$OUT/gbff_clean; mkdir -p "$CLEAN"
/home/shuaiw/miniconda3/envs/lovis4u/bin/python3 \
  /home/shuaiw/MODIFI/benchmark/kp_ECEs/clean_gbff.py "$GB" "$CLEAN"
GB=$CLEAN

echo "[pass1] feature/coord table"
$L4 -gb "$GB" -c "$CFG"  -alip -mmseqs-off -o "$OUT/lovis4u_pass1" --pdf-name pass1.pdf \
   > "$OUT/lovis4u_pass1.log" 2>&1 || { echo "pass1 failed"; tail -20 "$OUT/lovis4u_pass1.log"; exit 1; }

echo "[faf] map resistance genes"
python3 /home/shuaiw/MODIFI/benchmark/kp_ECEs/build_faf.py \
  "$OUT/lovis4u_pass1/feature_annotation_table.tsv" "$OUT/synteny_amrfinder.tsv" \
  "$OUT/synteny_faf.tsv" "$FIG/fig_infant_15_35_C_synteny_sourcedata.csv"
cp "$OUT/synteny_faf.tsv" "$FIG/synteny_faf.tsv"

echo "[pass2] final synteny figure"
$L4 -gb "$GB" -c "$CFG" -faf "$OUT/synteny_faf.tsv" -alip -hl -rol -oc \
   -o "$OUT/lovis4u_out" --pdf-name fig_infant_15_35_C_synteny.pdf \
   > "$OUT/lovis4u_pass2.log" 2>&1 || { echo "pass2 failed"; tail -25 "$OUT/lovis4u_pass2.log"; exit 1; }
cp "$OUT/lovis4u_out/fig_infant_15_35_C_synteny.pdf" "$FIG/"

echo "[png] rasterize"
if command -v pdftoppm >/dev/null 2>&1; then P=pdftoppm; else P="conda run -n mod pdftoppm"; fi
$P -png -r 150 "$FIG/fig_infant_15_35_C_synteny.pdf" "$FIG/fig_infant_15_35_C_synteny" >/dev/null 2>&1 \
  && ls "$FIG"/fig_infant_15_35_C_synteny*.png 2>/dev/null | head
echo "DONE"
