#!/bin/bash
#SBATCH --job-name=mix_10
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_10/modifi"

python /home/shuaiw/MODIFI/main.py \
    --work_dir "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_10/modifi/mix_10" \
    --unaligned_bam "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_10/mix_10.bam" \
    --whole_ref "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_10/mix_10.ref.fa" \
    --read_type hifi \
    --mge_file "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_10/mix_10.mge_list.csv" \
    --threads 64
