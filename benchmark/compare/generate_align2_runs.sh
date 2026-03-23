#!/usr/bin/env bash
# Emits one Slurm job script per reference under ./runs/ (same workload as batch_align2.sh).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_DIR="$SCRIPT_DIR/runs"
mkdir -p "$RUNS_DIR"

subreads_bam=/home/shuaiw/borg/XRSBK_20221007_S64018_PL100268287-1_C01.subreads.bam
ref_dir=/home/shuaiw/borg/paper/ipdsummary/subset_ref
ref_set=(test_50.fa test_100.fa test_200.fa test_300.fa test_400.fa test_500.fa)
outdir=/home/shuaiw/borg/paper/ipdsummary/soil_subset

for ref in "${ref_set[@]}"; do
    base="$(basename "$ref" .fa)"
    outfile="$RUNS_DIR/align2_${base}.sh"
    prefix="$outdir/$base"
    align_bam="${prefix}.align.subreads.bam"
    unsorted_bam="${prefix}.align.unsorted.bam"
    ref_path="$ref_dir/$ref"

    cat >"$outfile" <<EOF
#!/bin/bash
#SBATCH --job-name=align2_${base}
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "$outdir"

/usr/bin/time -v -o "${prefix}.pbmm2.align.time" \\
    ~/smrtlink/pbmm2 align "$ref_path" "$subreads_bam" "$unsorted_bam" \\
    --preset SUBREAD -j 64
samtools sort -T "$prefix" -@ 64 -o "$align_bam" "$unsorted_bam"
rm -f "$unsorted_bam"
samtools index "$align_bam"
~/smrtlink/pbindex "$align_bam"
EOF
    chmod +x "$outfile"
    echo "Wrote $outfile"
done
