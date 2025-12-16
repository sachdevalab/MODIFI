
                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_60/soil_60_methylation4 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_60
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_90/soil_90_methylation4 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_90
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_80/soil_80_methylation4 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_80
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_100/soil_100_methylation4 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2024.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_100
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_110/soil_110_methylation4 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2025.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_110
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_115/soil_115_methylation4 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_115
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_1/soil_1_methylation4 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_1
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/ece_LD/run/soil_2/soil_2_methylation4 \
                --unaligned_bam /home/shuaiw/borg/XRSBK_20221007_S64018_PL100268288-1_D01.ccs.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/ece_LD/22_ECE_and_Mp_seqs.fasta \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 3 \
                --min_ctg_cov 3 \
                --min_iden 0.97 \
                --min_frac 0.3 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
                --threads 64 " \
                --job-name=soil_2
            
