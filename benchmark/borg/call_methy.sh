python /home/shuaiw/MODIFI/main.py \
          --work_dir /home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_90 \
          --unaligned_bam  /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
          --whole_ref /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.fa \
          --read_type hifi \
            --min_len 1000 \
            --min_cov 3 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100  \
            --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
            --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
          --threads 10

python /home/shuaiw/MODIFI/main.py \
          --work_dir /home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_80 \
          --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
          --whole_ref /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_almost_complete_Black_borg_32_00.contigs.fa \
          --read_type hifi \
            --min_len 1000 \
            --min_cov 3 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100  \
            --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
            --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
          --threads 10

python /home/shuaiw/MODIFI/main.py \
          --work_dir /home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_60 \
          --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
          --whole_ref /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.contigs.fa \
          --read_type hifi \
            --min_len 1000 \
            --min_cov 3 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100  \
            --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
            --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
          --threads 10

python /home/shuaiw/MODIFI/main.py \
          --work_dir /home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_115 \
          --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
          --whole_ref /home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.contigs.fa \
          --read_type hifi \
            --min_len 1000 \
            --min_cov 3 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100  \
            --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
            --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
          --threads 10



python /home/shuaiw/MODIFI/main.py \
          --work_dir /home/shuaiw/borg/paper/borg_data/batch_export/soil_80_black \
          --unaligned_bam /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
          --whole_ref /home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.fasta \
          --read_type hifi \
            --min_len 1000 \
            --min_cov 3 \
            --min_frac 0.3 \
            --min_score 30 \
            --min_sites 100  \
            --kmer_mean_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.mean.dat \
            --kmer_num_db /home/shuaiw/borg/paper/run2/soil_1/soil_1_methylation4/control/control_db.up7.down3.num.dat \
          --threads 10