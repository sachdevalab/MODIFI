#!/usr/bin/env bash
set -euo pipefail

# Paths
MAIN="/home/shuaiw/Methy/main.py"
KMER_MEAN="/home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat"
KMER_NUM="/home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat"

# Samples
SAMPLES=(
  "NANO_2_INF1330004_4PB"
  "NANO_2_INF1340011_4PB"
)

# SAMPLES=(
#   "NANO_3_INF1330004_5PB"
#   "NANO_3_INF1330004_7PB"
#   "NANO_3_INF1330004_8PB"
#   "NANO_2_INF1340011_2PB"
#   "NANO_2_INF1340011_3PB"
#   "NANO_3_INF1340011_6PB"
#   "NANO_3_INF1340011_7PB"
# )

for sample in "${SAMPLES[@]}"; do
  bam="/home/shuaiw/borg/allison/batch/${sample}.align.ccs.bam"
  ref="/groups/banfield/projects/human/nano/1_assembly/pacbio_assembly/${sample}/${sample}_HR_HIFIASM_META_scaffold_min1000.fa"
  work="/home/shuaiw/borg/allison/batch2/${sample}_kmer"
  wrap=(python "$MAIN" --work_dir "$work" --whole_bam "$bam" --whole_ref "$ref" --read_type hifi --min_len 1000 --min_cov 5 --min_iden 0.97 --min_frac 0.3 --min_score 30 --min_sites 30 --kmer_mean_db "$KMER_MEAN" --kmer_num_db "$KMER_NUM" --threads 64)


  # Print commands instead of running them
  printf 'sbatch --partition standard --job-name=%s --wrap "%s"\n' "$sample" "${wrap[*]}"
  printf '\n'
done

