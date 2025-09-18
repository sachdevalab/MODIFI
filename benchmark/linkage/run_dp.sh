
            sbatch  --partition standard --wrap " /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p100.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p100 \
            --whole_bam /home/shuaiw/borg/paper/linkage/m64004_210929_143746.p100.bam \
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 --run_steps host" \
            --job-name=p100_pure
        

            sbatch  --partition standard --wrap " /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p10.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p10 \
            --whole_bam /home/shuaiw/borg/paper/linkage/m64004_210929_143746.p10.bam \
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 --run_steps host" \
            --job-name=p10_pure
        

            sbatch  --partition standard --wrap " /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p20.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p20 \
            --whole_bam /home/shuaiw/borg/paper/linkage/m64004_210929_143746.p20.bam \
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 --run_steps host" \
            --job-name=p20_pure
        

            sbatch  --partition standard --wrap " /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p30.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p30 \
            --whole_bam /home/shuaiw/borg/paper/linkage/m64004_210929_143746.p30.bam \
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 --run_steps host" \
            --job-name=p30_pure
        

            sbatch  --partition standard --wrap " /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p50.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p50 \
            --whole_bam /home/shuaiw/borg/paper/linkage/m64004_210929_143746.p50.bam \
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 --run_steps host" \
            --job-name=p50_pure
        

            sbatch  --partition standard --wrap " /usr/bin/time -v -o /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p05.time python /home/shuaiw/Methy/main.py \
            --work_dir /home/shuaiw/borg/paper/linkage/pure/m64004_210929_143746.p05 \
            --whole_bam /home/shuaiw/borg/paper/linkage/m64004_210929_143746.p05.bam \
            --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
            --read_type hifi \
            --min_len 1000 \
            --min_cov 1 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100 \
            --min_iden 0.97 \
            --mge_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list.tsv \
            --threads 64 --run_steps host" \
            --job-name=p05_pure
        
