#!/bin/bash
#SBATCH --job-name=mix_50
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_50/modifi"

python /home/shuaiw/MODIFI/main.py \
    --work_dir "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_50/modifi/mix_50" \
    --unaligned_bam "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_50/mix_50.bam" \
    --whole_ref "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_50/mix_50.ref.fa" \
    --read_type hifi \
    --mge_file "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_50/mix_50.mge_list.csv" \
    --threads 64
