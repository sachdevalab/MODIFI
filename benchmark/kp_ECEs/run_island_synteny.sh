#!/bin/bash
# Focused LoVis4u synteny of the two resistance islands of the infant_15_35_C cluster:
#   (1) copper/silver metal island (pco/sil +/- mer/ars/ter)   -> 11 members
#   (2) ars-integron AMR island (aadA1-sul1-qacEdelta1 / aph)  -> 4 members
# Each locus is windowed to its island (+/-4 kb); every resistance gene is labelled.
set -euo pipefail
OUT=/home/shuaiw/borg/revision/kp_eces/synteny
FIG=/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs
CFG=/home/shuaiw/MODIFI/benchmark/kp_ECEs/lovis4u_synteny.cfg
GB=$OUT/gbff_clean
L4="conda run -n lovis4u lovis4u"

# windows
python3 /home/shuaiw/MODIFI/benchmark/kp_ECEs/build_island_windows.py >/dev/null
# label EVERY resistance gene for the zoomed panels (room to show them all)
awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} {$4=1; print}' "$OUT/synteny_faf.tsv" > "$OUT/synteny_faf_all.tsv"

draw() {
  name=$1; winfile=$2
  $L4 -gb "$GB" -c "$CFG" -faf "$OUT/synteny_faf_all.tsv" -alip \
     -w $(cat "$OUT/$winfile") -hl -rol -oc -cl-owp \
     -o "$OUT/island_$name" --pdf-name "fig_infant_15_35_C_${name}_island.pdf" \
     > "$OUT/island_$name.log" 2>&1 || { echo "$name failed"; tail -20 "$OUT/island_$name.log"; exit 1; }
  cp "$OUT/island_$name/fig_infant_15_35_C_${name}_island.pdf" "$FIG/"
  if command -v pdftoppm >/dev/null 2>&1; then P=pdftoppm; else P="conda run -n mod pdftoppm"; fi
  $P -png -r 200 "$FIG/fig_infant_15_35_C_${name}_island.pdf" "$FIG/fig_infant_15_35_C_${name}_island" >/dev/null 2>&1 || true
  echo "wrote $name island figure"
}

draw metal metal_windows.txt
draw amr   amr_windows.txt
echo DONE
