
            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
                prefix=soil_60 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_60 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI.contigs.fa -j 64"  --job-name=soil_60
            

            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                prefix=soil_90 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_90 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI.contigs.fa -j 64"  --job-name=soil_90
            

            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
                prefix=soil_80 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_80 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI.contigs.fa -j 64"  --job-name=soil_80
            

            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2024.bam \
                prefix=soil_100 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_100 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI.contigs.fa -j 64"  --job-name=soil_100
            

            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2025.bam \
                prefix=soil_110 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_110 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_110cm_PACBIO-HIFI.contigs.fa -j 64"  --job-name=soil_110
            

            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
                prefix=soil_115 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_115 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI.contigs.fa -j 64"  --job-name=soil_115
            

            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam \
                prefix=soil_1 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_1 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.soil_1.fa -j 64"  --job-name=soil_1
            

            sbatch  --partition standard --wrap "snakemake -s soil_ggkbase_gtdb.smk --config \
                hifi_bam=/home/shuaiw/borg/XRSBK_20221007_S64018_PL100268288-1_D01.ccs.bam \
                prefix=soil_2 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_2 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_9_9_2021_34_2B_1_4m_PACBIO-HIFI_HIFIASM-META.soil_2.fa -j 64"  --job-name=soil_2
            
