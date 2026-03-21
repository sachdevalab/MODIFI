#!/usr/bin/env bash
# Emits one Slurm job script per reference under ./runs/ (same workload as batch_ipd2.sh).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS_DIR="$SCRIPT_DIR/runs"
mkdir -p "$RUNS_DIR"

ref_dir=/home/shuaiw/borg/contigs/
ref_set=(test_100.fa test_200.fa test_300.fa test_500.fa)
outdir=/home/shuaiw/borg/paper/ipdsummary/soil_1/
ipd_output=/home/shuaiw/borg/paper/ipdsummary/soil_1/ipd.out2/

for ref in "${ref_set[@]}"; do
    base="$(basename "$ref" .fa)"
    outfile="$RUNS_DIR/ipd2_${base}.sh"
    prefix="$outdir/$(basename "$ref" .fa)"
    ipd_prefix="$ipd_output/$(basename "$ref" .fa)"
    align_bam="$prefix.align.subreads.bam"
    ref_path="$ref_dir/$ref"

    cat >"$outfile" <<EOF
#!/bin/bash
#SBATCH --job-name=ipd2_${base}
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "$outdir"
mkdir -p "$ipd_output"

/usr/bin/time -v -o "$ipd_prefix.ipdSummary.time" \\
    ~/smrtlink/ipdSummary "$align_bam" \\
    --reference "$ref_path" \\
    --numWorkers 64 \\
    --gff "$ipd_prefix.gff" \\
    --csv "$ipd_prefix.csv" \\
    --methylFraction \\
    --outfile "$ipd_prefix.out"
EOF
    chmod +x "$outfile"
    echo "Wrote $outfile"
done
