#!/bin/bash
# run_e2e.sh — toy end-to-end de-novo run (Part H/C4 first test), terminal, 16 threads.
# assemble -> label contig origin (skani) -> group by origin -> CheckM2 -> geNomad -> MODIFI.
# Each step is guarded (skips if its output exists) so the driver is re-runnable.
set -uo pipefail
OUT=/home/shuaiw/borg/paper/simu_meta_dir/C4_toy
E2E=/home/shuaiw/MODIFI/benchmark/simu_meta/e2e
DB=/groups/diamond/databases/genomad/v1.7/
GGRENAME=/home/rohan/dev/pipeline/workflow/scripts/gg_rename_assembly.py
PY=/home/shuaiw/miniconda3/envs/modifi/bin/python
export PATH=/shared/software/bin:$PATH
T=16
cd "$OUT"

step(){ echo -e "\n=== [$(date +%H:%M)] $* ==="; }

# 1. de-novo co-assembly ----------------------------------------------------------
step "1. hifiasm_meta co-assembly"
if [ ! -f "$OUT/toy.hifiasm.p_ctg.gfa" ]; then
  hifiasm_meta --no-binning -o "$OUT/toy.hifiasm" -t $T "$OUT/toy.fq.gz"
fi
if [ ! -f "$OUT/toy.contigs.fa" ]; then
  awk '$1=="S"{printf ">%s\n%s\n",$2,$3}' "$OUT/toy.hifiasm.p_ctg.gfa" > "$OUT/toy.p_ctg.fa"
  if [ -r "$GGRENAME" ]; then
    $PY "$GGRENAME" -i "$OUT/toy.p_ctg.fa" -o "$OUT/toy.contigs.fa" -s toy
  else
    echo "  gg_rename not readable; falling back to plain rename (no _C/_L circular tag)"
    awk '/^>/{n++; print ">toy_"n"_L"; next} {print}' "$OUT/toy.p_ctg.fa" > "$OUT/toy.contigs.fa"
  fi
  samtools faidx "$OUT/toy.contigs.fa"
fi
echo "  contigs: $(grep -c '^>' "$OUT/toy.contigs.fa")"

# 2. contig origin ground truth (skani) -------------------------------------------
step "2. label contig origin (skani vs input genomes)"
[ -f "$OUT/contig_origin.tsv" ] || $PY "$E2E/label_contig_origin.py" \
    "$OUT/toy.contigs.fa" "$OUT/input_genomes" "$OUT/contig_origin.tsv" $T

# 3. group contigs by origin -> per-genome MAGs -> CheckM2 -------------------------
step "3. group by origin + CheckM2 per source genome"
if [ ! -f "$OUT/checkm2/quality_report.tsv" ]; then
  $PY - "$OUT/toy.contigs.fa" "$OUT/contig_origin.tsv" "$OUT/mags" <<'PYEOF'
import sys, os, pandas as pd
from collections import defaultdict
fa, lab, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)
origin = pd.read_csv(lab, sep="\t").set_index("contig")["origin_sample"].to_dict()
seqs, name = defaultdict(list), None
cur=[]
def flush():
    if name is not None:
        o = origin.get(name, "UNMAPPED")
        if o not in ("UNMAPPED",): seqs[o].append((name, "".join(cur)))
for line in open(fa):
    if line.startswith(">"):
        flush(); name=line[1:].strip().split()[0]; cur=[]
    else: cur.append(line.strip())
flush()
for o, recs in seqs.items():
    with open(os.path.join(outdir, f"{o}.fasta"),"w") as w:
        for n,s in recs: w.write(f">{n}\n{s}\n")
print(f"wrote {len(seqs)} per-origin MAG fastas")
PYEOF
  checkm2 predict --input "$OUT/mags" --output-directory "$OUT/checkm2" \
      --force -x .fasta --threads $T
fi
echo "  checkm2: $([ -f "$OUT/checkm2/quality_report.tsv" ] && echo OK || echo MISSING)"

# 4. geNomad ECE calling ----------------------------------------------------------
# geNomad ships as a Singularity container whose /opt/conda/bin (mmseqs+genomad) is NOT
# on PATH by default -> geNomad's which('mmseqs') fails. Bypass the wrapper and set PATH
# inside the container. /home and /groups are auto-bound (borg -> /groups/...).
step "4. geNomad end-to-end"
SIF=/shared/software/genomad/1.11.0/container/genomad.sif
if [ ! -d "$OUT/Genomad/toy.contigs_summary" ]; then
  /usr/local/bin/singularity exec --bind /home,/groups "$SIF" bash -c \
    "export PATH=/opt/conda/bin:\$PATH; genomad end-to-end --relaxed --enable-score-calibration --sensitivity 7.0 --threads $T --force-auto --cleanup '$OUT/toy.contigs.fa' '$OUT/Genomad' '$DB'"
fi
[ -f "$OUT/toy.mge.tsv" ] || $PY "$E2E/genomad_to_mge.py" \
    "$OUT/Genomad/toy.contigs_summary" "$OUT/toy.mge.tsv"

# 5. MODIFI (contig-level, NO --bin_file) -----------------------------------------
step "5. MODIFI host linkage (contig-level)"
if [ ! -f "$OUT/modifi/toy/host_summary.csv" ]; then
  mkdir -p "$OUT/modifi"
  $PY /home/shuaiw/MODIFI/main.py -o "$OUT/modifi/toy" -b "$OUT/toy.bam" \
      -r "$OUT/toy.contigs.fa" --read_type hifi --mge_file "$OUT/toy.mge.tsv" --threads $T
fi

# 6. score ------------------------------------------------------------------------
step "6. score: linkage success vs completeness"
$PY "$E2E/score_e2e.py" "$OUT"
echo -e "\n=== [$(date +%H:%M)] run_e2e.sh DONE ==="
