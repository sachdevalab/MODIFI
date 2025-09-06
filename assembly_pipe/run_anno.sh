
            #### number 4
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/XSDIA_20240410_R84050_PL8707-001_1-1-A01_bc2007-bc2007.hifi_reads.bam \
                prefix=cow_bioreactor_3 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_bioreactor_3 -j 64 --use-conda" \
                --job-name=cow_bioreactor_3 
            

            #### number 7
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow/RuReacBro_20230708_Comb_RF_LR.hifi_reads.bam \
                prefix=cow_1 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_1 -j 64 --use-conda" \
                --job-name=cow_1 
            

            #### number 11
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621429/ERR5621429.ccs.bam \
                prefix=ERR5621429_sludge \
                work_dir=/home/shuaiw/borg/paper/run2/ERR5621429_sludge -j 64 --use-conda" \
                --job-name=ERR5621429_sludge 
            

            #### number 12
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621430/ERR5621430.ccs.bam \
                prefix=ERR5621430_sludge \
                work_dir=/home/shuaiw/borg/paper/run2/ERR5621430_sludge -j 64 --use-conda" \
                --job-name=ERR5621430_sludge 
            

            #### number 17
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1310010_2G1_pacbio.bam \
                prefix=infant_5 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_5 -j 64 --use-conda" \
                --job-name=infant_5 
            

            #### number 23
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_3G1_pacbio.bam \
                prefix=infant_11 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_11 -j 64 --use-conda" \
                --job-name=infant_11 
            

            #### number 34
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340011_6G1_pacbio.bam \
                prefix=infant_22 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_22 -j 64 --use-conda" \
                --job-name=infant_22 
            

            #### number 37
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240045_6G1_pacbio.bam \
                prefix=infant_25 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_25 -j 64 --use-conda" \
                --job-name=infant_25 
            

            #### number 42
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
                prefix=soil_s3_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s3_1 -j 64 --use-conda" \
                --job-name=soil_s3_1 
            

            #### number 43
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                prefix=soil_s3_2 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s3_2 -j 64 --use-conda" \
                --job-name=soil_s3_2 
            

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
            
