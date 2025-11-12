
                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/run2//asthma_9/asthma_9_methylation3 \
                --whole_bam /home/shuaiw/borg/paper/run2//asthma_9/asthma_9.align.bam \
                --whole_ref /home/shuaiw/borg/paper/run2//asthma_9/asthma_9.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --run_steps control \
                --mge_file /home/shuaiw/borg/paper/run2//asthma_9/all_mge.tsv \
                --threads 64" \
                --job-name=sue_9
            
