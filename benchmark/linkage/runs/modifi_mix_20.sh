#!/bin/bash
#SBATCH --job-name=mix_20
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_20/modifi"

python /home/shuaiw/MODIFI/main.py \
    --work_dir "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_20/modifi/mix_20" \
    --unaligned_bam "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_20/mix_20.bam" \
    --whole_ref "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_20/mix_20.ref.fa" \
    --read_type hifi \
    --mge_file "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_20/mix_20.mge_list.csv" \
    --threads 64
