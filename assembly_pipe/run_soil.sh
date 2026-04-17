
sbatch --partition standard --job-name=soil_60 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_60/soil_60_methylation4 \
        -b /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"

sbatch --partition standard --job-name=soil_90 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_90/soil_90_methylation4 \
        -b /home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"

sbatch --partition standard --job-name=soil_80 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_80/soil_80_methylation4 \
        -b /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"

sbatch --partition standard --job-name=soil_100 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_100/soil_100_methylation4 \
        -b /home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2024.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"

sbatch --partition standard --job-name=soil_110 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_110/soil_110_methylation4 \
        -b /home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2025.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_110cm_PACBIO-HIFI.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"

sbatch --partition standard --job-name=soil_115 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_115/soil_115_methylation4 \
        -b /home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"

sbatch --partition standard --job-name=soil_1 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_1/soil_1_methylation4 \
        -b /home/shuaiw/borg/paper/aws/soil/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"

sbatch --partition standard --job-name=soil_2 --wrap "python /home/shuaiw/MODIFI/main.py \
        -o /home/shuaiw/borg/paper/gg_run4/soil_2/soil_2_methylation4 \
        -b /home/shuaiw/borg/XRSBK_20221007_S64018_PL100268288-1_D01.ccs.bam \
        -r /home/shuaiw/borg/paper/curated_genome/unique/SR-VP_9_9_2021_34_2B_1_4m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
        --read_type hifi --min_len 1000 --min_cov 3 --min_frac 0.3 --min_score 30 \
        --min_sites 100 --min_ctg_cov 2 \
        --mge_file /home/shuaiw/MODIFI/benchmark/borg/klingon/klingon_contigs.txt \
        --threads 64 --no-clean"
