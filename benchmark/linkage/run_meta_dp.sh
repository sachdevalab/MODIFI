#!/bin/bash
#SBATCH --job-name=plex_soil 
 #SBATCH --partition=standard

            /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p10.run.time python /home/shuaiw/Methy/main.py             --work_dir /home/shuaiw/borg/paper/linkage/meta/             --whole_bam /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p10.align.bam             --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa             --read_type hifi             --min_len 1000             --max_NM 10             --min_cov 1             --min_frac 0.4             --min_score 30             --min_sites 30             --clean             --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list
        

            /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p20.run.time python /home/shuaiw/Methy/main.py             --work_dir /home/shuaiw/borg/paper/linkage/meta/             --whole_bam /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p20.align.bam             --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa             --read_type hifi             --min_len 1000             --max_NM 10             --min_cov 1             --min_frac 0.4             --min_score 30             --min_sites 30             --clean             --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list
        

            /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p30.run.time python /home/shuaiw/Methy/main.py             --work_dir /home/shuaiw/borg/paper/linkage/meta/             --whole_bam /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p30.align.bam             --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa             --read_type hifi             --min_len 1000             --max_NM 10             --min_cov 1             --min_frac 0.4             --min_score 30             --min_sites 30             --clean             --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list
        

            /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p50.run.time python /home/shuaiw/Methy/main.py             --work_dir /home/shuaiw/borg/paper/linkage/meta/             --whole_bam /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p50.align.bam             --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa             --read_type hifi             --min_len 1000             --max_NM 10             --min_cov 1             --min_frac 0.4             --min_score 30             --min_sites 30             --clean             --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list
        

            /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p05.run.time python /home/shuaiw/Methy/main.py             --work_dir /home/shuaiw/borg/paper/linkage/meta/             --whole_bam /home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.p05.align.bam             --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa             --read_type hifi             --min_len 1000             --max_NM 10             --min_cov 1             --min_frac 0.4             --min_score 30             --min_sites 30             --clean             --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list
        
