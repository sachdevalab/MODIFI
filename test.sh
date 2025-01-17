sbatch --partition memory --wrap "/usr/bin/time -v -o slurm.human.time python get_baseline3.py"
