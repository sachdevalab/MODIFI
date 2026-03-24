#!/bin/bash
#SBATCH --job-name=mix_30
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_30/modifi"

python /home/shuaiw/MODIFI/main.py \
    --work_dir "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_30/modifi/mix_30" \
    --unaligned_bam "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_30/mix_30.bam" \
    --whole_ref "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_30/mix_30.ref.fa" \
    --read_type hifi \
    --mge_file "/home/shuaiw/borg/paper/linkage/mixed_isolates/mix_30/mix_30.mge_list.csv" \
    --threads 64
