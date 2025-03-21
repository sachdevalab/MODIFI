sbatch  --partition standard --wrap "metabat --inFile /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged.fa\
        --outFile /home/shuaiw/borg/maxbat\
        -l --seed 8 -s 2000" --job-name=maxbat
        
