#!/bin/bash
# Build C2 strain-depth levels K=2,3,4,all (nested supersets of bg_300) in the terminal
# using the parallel builder (64 threads), then submit each MODIFI run to SLURM.
# Requires bg_300 already built (build_c2_strain.py reads its manifest as the base).
cd /home/shuaiw/MODIFI/benchmark/simu_meta
ENVBIN=/home/shuaiw/miniconda3/envs/modifi/bin; MAIN=/home/shuaiw/MODIFI/main.py
SM=/home/shuaiw/MODIFI/benchmark/simu_meta; DIRROOT=/home/shuaiw/borg/paper/simu_meta_dir/C1
export PATH="$ENVBIN:$PATH"
LOG=$DIRROOT/build_c2.log; : > "$LOG"

for K in 2 3 4 all; do
  label=bg300_k$K; d=$DIRROOT/$label
  echo "$(date) building $label ..." >> "$LOG"
  python3 build_c2_strain.py --K "$K" --threads 64 >> "$LOG" 2>&1
  if [ -f "$d/$label.bam" ]; then
    script=$SM/cmds/job_$label.sh
    cat > "$script" <<EOF
#!/bin/bash
#SBATCH --job-name=sm_$label
#SBATCH --partition=standard
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
set -euo pipefail
export PATH="$ENVBIN:\$PATH"
mkdir -p "$d/modifi"
$ENVBIN/python $MAIN -o "$d/modifi/$label" -b "$d/$label.bam" -r "$d/$label.ref.fa" \\
    --read_type hifi --mge_file "$d/$label.mge_list.tsv" --threads 32
EOF
    chmod +x "$script"; job=$(sbatch --parsable "$script")
    echo "$(date) submitted $label MODIFI -> job $job" >> "$LOG"
  else
    echo "$(date) BUILD FAILED for $label" >> "$LOG"
  fi
done
echo "$(date) === C2 build+submit done ===" >> "$LOG"
