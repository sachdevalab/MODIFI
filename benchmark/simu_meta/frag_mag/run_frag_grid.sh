#!/bin/bash
# run_frag_grid.sh — Step 2 of the synthetic Part H: completeness x contamination grid.
# PRE-REQ: Step 1 (MODIFI on fragmented ladder_58) has produced per-fragment profiles in
# $WORK. For each (completeness, contamination) cell we build a --bin_file (each host MAG =
# a subset of its fragments + foreign fragments) and run LINKAGE ONLY (--run_steps host),
# reusing the cached per-fragment methylation. Cheap: no re-alignment / re-calling.
set -uo pipefail
FM=/home/shuaiw/MODIFI/benchmark/simu_meta/frag_mag
D=/home/shuaiw/borg/paper/simu_meta_dir/C1/ladder58_frag
SRC=/home/shuaiw/borg/paper/simu_meta_dir/C1/ladder_58_rep2
WORK=$D/modifi/ladder58_frag                    # Step 1 output (has profiles/, motif_profile.csv)
PY=/home/shuaiw/miniconda3/envs/modifi/bin/python
MAIN=/home/shuaiw/MODIFI/main.py
T=${1:-32}
mkdir -p "$D/bins" "$D/results"

if [ ! -f "$WORK/motif_profile.csv" ]; then
  echo "ERROR: Step 1 profiles not found in $WORK — wait for the fragmentation MODIFI run."; exit 1
fi

for comp in 50 60 70 80 90 100; do
  for contam in 0 5 10 15; do
    tag="c${comp}_x${contam}"
    res="$D/results/${tag}.host_summary.csv"
    [ -f "$res" ] && { echo "skip $tag (done)"; continue; }
    binf="$D/bins/${tag}.bin"
    $PY "$FM/make_bins.py" "$D/fragment_map.tsv" "$binf" "$comp" "$contam"
    # MODIFI writes correct per-ECE bin predictions (hosts/) then crashes on its final
    # host_summary write in this bin setup -> tolerate the non-zero exit, aggregate hosts/ ourselves.
    $PY "$MAIN" --run_steps host -o "$WORK" -b "$SRC/ladder_58_rep2.bam" \
        -r "$D/ladder58_frag.ref.fa" --mge_file "$SRC/ladder_58_rep2.mge_list.tsv" \
        --bin_file "$binf" --threads "$T" > "$D/results/${tag}.linkage.log" 2>&1 || true
    $PY "$FM/agg_hosts.py" "$WORK/hosts" "$res" && echo "  $tag done ($(($(wc -l <"$res")-1)) ECEs)"
  done
done
$PY "$FM/score_mag.py" "$D/results"
echo "=== frag grid DONE ==="
