sbatch  --partition standard --wrap 'bash batch/run_isolation_part_0.sh' --job-name=iso_0
sbatch  --partition standard --wrap 'bash batch/run_isolation_part_1.sh' --job-name=iso_1
sbatch  --partition standard --wrap 'bash batch/run_isolation_part_2.sh' --job-name=iso_2
