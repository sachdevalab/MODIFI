

            snakemake -s soil_ggkbase_host.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2024.bam \
                prefix=soil_100 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_100 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI.contigs.fa -j 64
            

            snakemake -s soil_ggkbase_host.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
                prefix=soil_115 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_115 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI.contigs.fa -j 64


            snakemake -s soil_ggkbase_host.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                prefix=soil_90 \
                work_dir=/home/shuaiw/borg/paper/gg_run2/soil_90 \
                ref=/home/shuaiw/borg/paper/curated_genome/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI.contigs.fa -j 64
            
