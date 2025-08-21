

            #### number 43
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                prefix=soil_s3_2 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s3_2 -j 64 --use-conda" \
                --job-name=soil_s3_2 
            

            #### number 44
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
                prefix=soil_s4_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s4_1 -j 64 --use-conda" \
                --job-name=soil_s4_1 
            

            #### number 46
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2025.bam \
                prefix=soil_s1_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s1_1 -j 64 --use-conda" \
                --job-name=soil_s1_1 
            

            #### number 47
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
                prefix=soil_s1_2 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s1_2 -j 64 --use-conda" \
                --job-name=soil_s1_2 
            

            #### number 48
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam \
                prefix=soil_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_1 -j 64 --use-conda" \
                --job-name=soil_1 
            

            #### number 49
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR14074352/SRR14074352.ccs.bam \
                prefix=SRR14074352_human \
                work_dir=/home/shuaiw/borg/paper/run2/SRR14074352_human -j 64 --use-conda" \
                --job-name=SRR14074352_human 
            
