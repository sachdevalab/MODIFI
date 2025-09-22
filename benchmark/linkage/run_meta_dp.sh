#!/bin/bash
#SBATCH --job-name=plex_infant 
 #SBATCH --partition=standard

            sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.run.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100/ \
            --whole_bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p100.align.bam \
            --whole_ref /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 " \
            --job-name=p100
        

            sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.run.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10/ \
            --whole_bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p10.align.bam \
            --whole_ref /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 " \
            --job-name=p10
        

            sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.run.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20/ \
            --whole_bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p20.align.bam \
            --whole_ref /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 " \
            --job-name=p20
        

            sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.run.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30/ \
            --whole_bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p30.align.bam \
            --whole_ref /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 " \
            --job-name=p30
        

            sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.run.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50/ \
            --whole_bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p50.align.bam \
            --whole_ref /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 " \
            --job-name=p50
        

            sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.run.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05/ \
            --whole_bam /home/shuaiw/borg/paper/linkage/meta_infant_14/m64004_210929_143746.p05.align.bam \
            --whole_ref /home/shuaiw/borg/paper/linkage/meta_infant_14/96plex_infant.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.4 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 " \
            --job-name=p05
        
