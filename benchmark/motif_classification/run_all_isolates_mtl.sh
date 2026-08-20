#!/bin/bash
# Batch: run Nanomotif MTase-linker on every isolate that has BOTH an assembly and detected
# motifs. Terminal run, capped at 64 threads total (32 parallel jobs x 2 threads each).
# Idempotent: isolates with an existing mtase_assignment_table.tsv are skipped.
set -uo pipefail

source /home/shuaiw/miniconda3/etc/profile.d/conda.sh 2>/dev/null
export PATH=/home/shuaiw/miniconda3/envs/nanomotif/bin:$PATH

BATCH=/home/shuaiw/borg/paper/isolation/batch2_results
RUNNER=/home/shuaiw/MODIFI/benchmark/motif_classification/run_isolate_mtl.sh
LIST=/home/shuaiw/borg/revision/motif_class/isolate_mtl/isolate_list.txt
mkdir -p /home/shuaiw/borg/revision/motif_class/isolate_mtl

# Build the ACC list: assembly AND all.motifs.csv both present.
comm -12 \
  <(find "$BATCH" -maxdepth 2 -name "*.hifiasm.p_ctg.rename.fa" 2>/dev/null \
      | sed -E 's#.*/([^/]+)/[^/]+$#\1#' | sort -u) \
  <(find "$BATCH" -maxdepth 3 -name "all.motifs.csv" 2>/dev/null \
      | sed -E 's#.*/([^/]+)/[^/]+_methylation4/all.motifs.csv#\1#' | sort -u) \
  > "$LIST"
echo "isolates to process: $(wc -l < "$LIST")"

# 32 concurrent jobs x 2 threads = 64 threads.
cat "$LIST" | xargs -P 32 -I{} bash "$RUNNER" {} 2 >/dev/null 2>&1

echo "ALL DONE. tables: $(find /home/shuaiw/borg/revision/motif_class/isolate_mtl -name mtase_assignment_table.tsv | wc -l)"
