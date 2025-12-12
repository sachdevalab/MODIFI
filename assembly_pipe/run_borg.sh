
                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_s3_1/soil_s3_1_methylation3 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_42
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_s3_2/soil_s3_2_methylation3 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_43
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_s4_1/soil_s4_1_methylation3 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_44
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_s4_2/soil_s4_2_methylation3 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2024.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_45
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_s1_1/soil_s1_1_methylation3 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2025.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_46
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_s1_2/soil_s1_2_methylation3 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_47
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_1/soil_1_methylation3 \
                --unaligned_bam /home/shuaiw/borg/paper/aws/soil/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_48
            

                sbatch  --partition standard --wrap "python /home/shuaiw/mGlu/main.py \
                --work_dir /home/shuaiw/borg/paper/borg_data/borg_for2/soil_2/soil_2_methylation3 \
                --unaligned_bam /home/shuaiw/borg/XRSBK_20221007_S64018_PL100268288-1_D01.ccs.bam \
                --whole_ref /home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa \
                --read_type hifi \
                --min_len 1000 \
                --min_cov 2 \
                --min_ctg_cov 2 \
                --min_iden 0.95 \
                --min_frac 0.1 \
                --min_score 30 \
                --min_sites 100 \
                --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.mean.dat \
                --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation3/control/control_db.up7.down3.num.dat \
                --mge_file /home/shuaiw/borg/paper/borg_data/align/borg.tsv \
                --threads 64 --visu_ipd --detect_misassembly" \
                --job-name=borg_49
            
