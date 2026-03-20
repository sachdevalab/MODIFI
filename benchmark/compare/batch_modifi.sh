#!/bin/bash
#SBATCH --job-name=modifi_28        # Job name
#SBATCH --partition=standard        # Partition name

set -euo pipefail

ref_dir=/home/shuaiw/borg/contigs/
ref_set=(test_100.fa test_200.fa test_300.fa test_500.fa)
modifi_out_root=/home/shuaiw/borg/paper/ipdsummary/soil_1/modifi.out/
threads="${SLURM_CPUS_ON_NODE:-64}"

mkdir -p "$modifi_out_root"

for ref in "${ref_set[@]}"; do
    base="$(basename "$ref" .fa)"
    ref_path="$ref_dir/$ref"
    work_dir="$modifi_out_root/$base"

    align_bam="/home/shuaiw/borg/paper/ipdsummary/soil_1/$base.align.subreads.bam"
    

    mkdir -p "$work_dir"

    /usr/bin/time -v -o "$work_dir/modifi.host.time" \
        python /home/shuaiw/MODIFI/main.py \
        --work_dir "$work_dir" \
        --whole_bam "$align_bam" \
        --whole_ref "$ref_path" \
        --read_type subreads \
        --run_steps split load control compare \
        --threads $threads
    # break
done
