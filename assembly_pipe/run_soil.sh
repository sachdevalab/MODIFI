
                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_60/soil_60_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_60

            

                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_90/soil_90_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_90

            

                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_80/soil_80_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_80

            

                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_100/soil_100_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_100

            

                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_110/soil_110_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_110cm_PACBIO-HIFI.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_110

            

                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_115/soil_115_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_115

            

                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_1/soil_1_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_1

            

                sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
                        -o /home/shuaiw/borg/paper/gg_run3/soil_2/soil_2_methylation4 \
                        --aligned_bam xxx \
                        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_9_9_2021_34_2B_1_4m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
                        --read_type hifi \
                        --min_len 1000 \
                        --min_cov 3 \
                        --min_frac 0.3 \
                        --min_score 30 \
                        --min_sites 100  \
                        --mge_file /home/shuaiw/borg/paper/natasha/klingon.genome.list \
                        --threads 64 --run_steps host --no-clean" \
                        --job-name=soil_2

            
