 sbatch  --partition standard --wrap "/usr/bin/time -v -o ipdSummary.time snakemake -j 64"  --job-name=ipdSum 
