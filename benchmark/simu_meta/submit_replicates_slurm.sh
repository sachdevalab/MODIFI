#!/bin/bash
# Generate + submit one self-contained SLURM job per replicate community.
# Each job BUILDS the merged community (samtools) AND runs MODIFI on the compute node,
# so nothing heavy runs on the login node. Replicate 1 = existing seed-42 runs;
# this handles replicates 2-5 (seeds 43-46). Throttled to <=9 concurrent jobs.
# The submitter loop is lightweight (sbatch/squeue/sleep) and runs detached.
set -u
SM=/home/shuaiw/MODIFI/benchmark/simu_meta
ENVBIN=/home/shuaiw/miniconda3/envs/modifi/bin
MAIN=/home/shuaiw/MODIFI/main.py
DIRROOT=/home/shuaiw/borg/paper/simu_meta_dir/C1
GEN=$SM/cmds
LOG=$DIRROOT/reps_slurm.log
mkdir -p "$GEN"; : > "$LOG"

njobs() { squeue -u shuaiw -h 2>/dev/null | wc -l; }
already_queued() { squeue -u shuaiw -h -o '%j' 2>/dev/null | grep -qx "sm_$1"; }

# write a combined build+run SLURM script; $1=label  $2=build command (relative to $SM)
# Submit ALL jobs at once; SLURM runs <=10 and queues the rest. Skip if already queued
# or already finished (host_summary present).
emit_and_submit() {
  local label="$1"; local buildcmd="$2"
  local dir="$DIRROOT/$label"
  local script="$GEN/job_$label.sh"
  if already_queued "$label"; then echo "$(date) skip $label (already queued)" >> "$LOG"; return; fi
  if [ -f "$dir/modifi/$label/host_summary.csv" ]; then echo "$(date) skip $label (done)" >> "$LOG"; return; fi
  cat > "$script" <<EOF
#!/bin/bash
#SBATCH --job-name=sm_$label
#SBATCH --partition=standard
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
set -euo pipefail
export PATH="$ENVBIN:\$PATH"
cd $SM
# --- build merged community if not already present ---
if [ ! -f "$dir/$label.bam" ]; then
  $buildcmd
fi
mkdir -p "$dir/modifi"
# --- run MODIFI ---
$ENVBIN/python $MAIN \\
    -o "$dir/modifi/$label" \\
    -b "$dir/$label.bam" \\
    -r "$dir/$label.ref.fa" \\
    --read_type hifi \\
    --mge_file "$dir/$label.mge_list.tsv" \\
    --threads 32
EOF
  chmod +x "$script"
  local job; job=$(sbatch --parsable "$script")
  echo "$(date) submitted $label -> job $job (jobs now $(njobs))" >> "$LOG"
}

echo "$(date) === replicate SLURM submitter start ===" >> "$LOG"

# --- ladder replicates 2-5 (each size builds itself from the seed's nested base via --only) ---
for rep in 2 3 4 5; do
  seed=$((41+rep))
  for s in 10 25 40 58; do
    emit_and_submit "ladder_${s}_rep$rep" \
      "python3 build_community.py ladder --sizes 10,25,40,58 --seed $seed --tag rep$rep --only $s --threads 32"
  done
done

# --- background communities replicates 2-5 ---
emit_bg() {  # $1=base $2=nsp $3=k $4=nbg
  for rep in 2 3 4 5; do
    seed=$((41+rep))
    emit_and_submit "${1}_rep$rep" \
      "python3 build_community.py community --n-species $2 --strains-per-species $3 --n-background $4 --label $1 --tag rep$rep --seed $seed --threads 32"
  done
}
# C1 = complexity only: donors are ALWAYS 1 strain/species (strains-per-species=1);
# background provides the extra genomes (background may repeat species/strains).
emit_bg bg_80  40 1 40
emit_bg bg_150 58 1 95
emit_bg bg_300 58 1 246

echo "$(date) === all replicate jobs submitted ===" >> "$LOG"
