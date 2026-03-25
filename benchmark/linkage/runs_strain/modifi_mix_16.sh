#!/bin/bash
#SBATCH --job-name=mix_16
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/linkage/mixed_isolates_strain/mix_16/modifi"

python /home/shuaiw/MODIFI/main.py \
    --work_dir "/home/shuaiw/borg/paper/linkage/mixed_isolates_strain/mix_16/modifi/mix_16" \
    --unaligned_bam "/home/shuaiw/borg/paper/linkage/mixed_isolates_strain/mix_16/mix_16.bam" \
    --whole_ref "/home/shuaiw/borg/paper/linkage/mixed_isolates_strain/mix_16/mix_16.ref.fa" \
    --read_type hifi \
    --mge_file "/home/shuaiw/borg/paper/linkage/mixed_isolates_strain/mix_16/mix_16.mge_list.csv" \
    --threads 64
