 sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p05.run.time python /home/shuaiw/Methy/main.py 
\
            --work_dir /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p05/ \
            --whole_bam /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p05.align.bam \
            --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa \
            --read_type hifi \
            --min_len 1000 \
            --max_NM 2000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 30 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 "\                                                                                                      
            --job-name=p05
