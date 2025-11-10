
                # sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                # --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_1/asthma_1_methylation3 \
                # --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_1/asthma_1.align.bam \
                # --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_1/asthma_1.hifiasm.p_ctg.rename.fa \
                # --read_type hifi \
                # --min_len 2000 \
                # --min_cov 3 \
                # --min_iden 0.97 \
                # --min_frac 0.3 \
                # --min_score 30 \
                # --min_sites 100 \
                # --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_1/all_mge.tsv \
                # --threads 64" \
                # --job-name=sue_1
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_2/asthma_2_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_2/asthma_2.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_2/asthma_2.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_2/all_mge.tsv \
                --threads 64" \
                --job-name=sue_2
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_3/asthma_3_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_3/asthma_3.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_3/asthma_3.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_3/all_mge.tsv \
                --threads 64" \
                --job-name=sue_3
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_4/asthma_4_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_4/asthma_4.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_4/asthma_4.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_4/all_mge.tsv \
                --threads 64" \
                --job-name=sue_4
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_5/asthma_5_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_5/asthma_5.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_5/asthma_5.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_5/all_mge.tsv \
                --threads 64" \
                --job-name=sue_5
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_6/asthma_6_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_6/asthma_6.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_6/asthma_6.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_6/all_mge.tsv \
                --threads 64" \
                --job-name=sue_6
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_7/asthma_7_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_7/asthma_7.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_7/asthma_7.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_7/all_mge.tsv \
                --threads 64" \
                --job-name=sue_7
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_8/asthma_8_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_8/asthma_8.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_8/asthma_8.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_8/all_mge.tsv \
                --threads 64" \
                --job-name=sue_8
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_9/asthma_9_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_9/asthma_9.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_9/asthma_9.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_9/all_mge.tsv \
                --threads 64" \
                --job-name=sue_9
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_10/asthma_10_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_10/asthma_10.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_10/asthma_10.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_10/all_mge.tsv \
                --threads 64" \
                --job-name=sue_10
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_11/asthma_11_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_11/asthma_11.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_11/asthma_11.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_11/all_mge.tsv \
                --threads 64" \
                --job-name=sue_11
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_12/asthma_12_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_12/asthma_12.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_12/asthma_12.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_12/all_mge.tsv \
                --threads 64" \
                --job-name=sue_12
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_13/asthma_13_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_13/asthma_13.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_13/asthma_13.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_13/all_mge.tsv \
                --threads 64" \
                --job-name=sue_13
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_14/asthma_14_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_14/asthma_14.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_14/asthma_14.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_14/all_mge.tsv \
                --threads 64" \
                --job-name=sue_14
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_15/asthma_15_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_15/asthma_15.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_15/asthma_15.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_15/all_mge.tsv \
                --threads 64" \
                --job-name=sue_15
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_16/asthma_16_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_16/asthma_16.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_16/asthma_16.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_16/all_mge.tsv \
                --threads 64" \
                --job-name=sue_16
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_17/asthma_17_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_17/asthma_17.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_17/asthma_17.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_17/all_mge.tsv \
                --threads 64" \
                --job-name=sue_17
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_18/asthma_18_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_18/asthma_18.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_18/asthma_18.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_18/all_mge.tsv \
                --threads 64" \
                --job-name=sue_18
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_19/asthma_19_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_19/asthma_19.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_19/asthma_19.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_19/all_mge.tsv \
                --threads 64" \
                --job-name=sue_19
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_20/asthma_20_methylation3 \
                --whole_bam /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_20/asthma_20.align.bam \
                --whole_ref /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_20/asthma_20.hifiasm.p_ctg.rename.fa \
                --read_type hifi \
                --min_len 2000 \
                --min_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --mge_file /groups/banfield/projects/multienv/methylation_temp/asthma_run//asthma_20/all_mge.tsv \
                --threads 64" \
                --job-name=sue_20
            
