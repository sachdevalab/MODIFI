#!/bin/bash
#SBATCH -p standard
#SBATCH -J ocra_bax2bam
#SBATCH --qos unlimit-submit-20-run
#SBATCH -o /home/shuaiw/borg/revision/ocra_5mC/ocra_bax2bam.slurm.out
set -euo pipefail
source /home/shuaiw/miniconda3/etc/profile.d/conda.sh
conda activate bax2bam
export PATH=/home/shuaiw/miniconda3/envs/bax2bam/bin:$PATH
SM=/home/shuaiw/smrtlink
R=/home/shuaiw/borg/revision/ocra_5mC
cd "$R"
M34=m170105_145115_42200_c101146242550000001823274307071741
M35=m170105_102734_42200_c101146242550000001823274307071740

# Convert each SMRT cell's 3 bax.h5 -> <SRR>.subreads.bam (keeps ip/pw kinetics)
bax2bam -o SRR8580034 \
  SRR8580034/${M34}_s1_p0.1.bax.h5 SRR8580034/${M34}_s1_p0.2.bax.h5 SRR8580034/${M34}_s1_p0.3.bax.h5
bax2bam -o SRR8580035 \
  SRR8580035/${M35}_s1_p0.1.bax.h5 SRR8580035/${M35}_s1_p0.2.bax.h5 SRR8580035/${M35}_s1_p0.3.bax.h5

# Merge the two cells for the deep (main) run
$SM/pbmerge -o ocra.subreads.bam SRR8580034.subreads.bam SRR8580035.subreads.bam
$SM/pbindex ocra.subreads.bam

echo "=== kinetics check (should show ip:B / pw:B) ==="
$SM/samtools view SRR8580034.subreads.bam 2>/dev/null | head -1 | grep -o -E 'ip:B:[A-Za-z]|pw:B:[A-Za-z]' | sort -u
echo "=== read counts ==="
for b in SRR8580034 SRR8580035 ocra; do echo -n "$b.subreads.bam: "; $SM/samtools view -c $b.subreads.bam; done
echo "BAX2BAM+MERGE DONE"
