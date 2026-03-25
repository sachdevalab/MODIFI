#!/bin/bash
#SBATCH --job-name=mix_40
#SBATCH --partition=standard

set -euo pipefail

mkdir -p "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_40/modifi"

python /home/shuaiw/MODIFI/main.py \
    --work_dir "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_40/modifi/mix_40" \
    --unaligned_bam "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_40/mix_40.bam" \
    --whole_ref "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_40/mix_40.ref.fa" \
    --read_type hifi \
    --mge_file "/home/shuaiw/borg/paper/linkage/mixed_isolates2/mix_40/mix_40.mge_list.csv" \
    --threads 64
