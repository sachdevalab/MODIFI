
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
            

            #### number 10
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621427/ERR5621427.ccs.bam \
                prefix=ERR5621427_sludge \
                work_dir=/home/shuaiw/borg/paper/run2/ERR5621427_sludge -j 64 --use-conda" \
                --job-name=ERR5621427_sludge 
            

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
                hifi_bam=/home/shuaiw/borg/XRSBK_20221007_S64018_PL100268288-1_D01.ccs.bam \
                prefix=soil_2 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_2 -j 64 --use-conda" \
                --job-name=soil_2 
            

            #### number 50
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR14074352/SRR14074352.ccs.bam \
                prefix=SRR14074352_human \
                work_dir=/home/shuaiw/borg/paper/run2/SRR14074352_human -j 64 --use-conda" \
                --job-name=SRR14074352_human 
            

            #### number 53
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208058140_12Mo_LR/PACBIO_DATA/TIPS_208058140_12Mo_LR.hifi_reads.bam \
                prefix=asthma_1 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_1 -j 64 --use-conda" \
                --job-name=asthma_1 
            

            #### number 54
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208058140_1Mo_LR/PACBIO_DATA/TIPS_208058140_1Mo_LR.hifi_reads.bam \
                prefix=asthma_2 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_2 -j 64 --use-conda" \
                --job-name=asthma_2 
            

            #### number 55
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208062141_12Mo_LR/PACBIO_DATA/TIPS_208062141_12Mo_LR.hifi_reads.bam \
                prefix=asthma_3 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_3 -j 64 --use-conda" \
                --job-name=asthma_3 
            

            #### number 56
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208062141_1Mo_LR/PACBIO_DATA/TIPS_208062141_1Mo_LR.hifi_reads.bam \
                prefix=asthma_4 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_4 -j 64 --use-conda" \
                --job-name=asthma_4 
            

            #### number 57
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208075159_12Mo_LR/PACBIO_DATA/TIPS_208075159_12Mo_LR.hifi_reads.bam \
                prefix=asthma_5 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_5 -j 64 --use-conda" \
                --job-name=asthma_5 
            

            #### number 58
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208075159_1Mo_LR/PACBIO_DATA/TIPS_208075159_1Mo_LR.hifi_reads.bam \
                prefix=asthma_6 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_6 -j 64 --use-conda" \
                --job-name=asthma_6 
            

            #### number 59
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208107157_12Mo_LR/PACBIO_DATA/TIPS_208107157_12Mo_LR.hifi_reads.bam \
                prefix=asthma_7 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_7 -j 64 --use-conda" \
                --job-name=asthma_7 
            

            #### number 60
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208107157_1Mo_LR/PACBIO_DATA/TIPS_208107157_1Mo_LR.hifi_reads.bam \
                prefix=asthma_8 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_8 -j 64 --use-conda" \
                --job-name=asthma_8 
            

            #### number 61
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209129175_12Mo_LR/PACBIO_DATA/TIPS_209129175_12Mo_LR.hifi_reads.bam \
                prefix=asthma_9 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_9 -j 64 --use-conda" \
                --job-name=asthma_9 
            

            #### number 62
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209129175_1Mo_LR/PACBIO_DATA/TIPS_209129175_1Mo_LR.hifi_reads.bam \
                prefix=asthma_10 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_10 -j 64 --use-conda" \
                --job-name=asthma_10 
            

            #### number 63
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209134177_12Mo_LR/PACBIO_DATA/TIPS_209134177_12Mo_LR.hifi_reads.bam \
                prefix=asthma_11 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_11 -j 64 --use-conda" \
                --job-name=asthma_11 
            

            #### number 64
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209134177_1Mo_LR/PACBIO_DATA/TIPS_209134177_1Mo_LR.hifi_reads.bam \
                prefix=asthma_12 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_12 -j 64 --use-conda" \
                --job-name=asthma_12 
            

            #### number 65
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209158193_12Mo_LR/PACBIO_DATA/TIPS_209158193_12Mo_LR.hifi_reads.bam \
                prefix=asthma_13 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_13 -j 64 --use-conda" \
                --job-name=asthma_13 
            

            #### number 66
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209158193_1Mo_LR/PACBIO_DATA/TIPS_209158193_1Mo_LR.hifi_reads.bam \
                prefix=asthma_14 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_14 -j 64 --use-conda" \
                --job-name=asthma_14 
            

            #### number 67
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_210173216_12Mo_LR/PACBIO_DATA/TIPS_210173216_12Mo_LR.hifi_reads.bam \
                prefix=asthma_15 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_15 -j 64 --use-conda" \
                --job-name=asthma_15 
            

            #### number 68
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_210173216_1Mo_LR/PACBIO_DATA/TIPS_210173216_1Mo_LR.hifi_reads.bam \
                prefix=asthma_16 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_16 -j 64 --use-conda" \
                --job-name=asthma_16 
            

            #### number 69
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_211245247_12Mo_LR/PACBIO_DATA/TIPS_211245247_12Mo_LR.hifi_reads.bam \
                prefix=asthma_17 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_17 -j 64 --use-conda" \
                --job-name=asthma_17 
            

            #### number 70
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_211245247_1Mo_LR/PACBIO_DATA/TIPS_211245247_1Mo_LR.hifi_reads.bam \
                prefix=asthma_18 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_18 -j 64 --use-conda" \
                --job-name=asthma_18 
            

            #### number 71
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_212256257_12Mo_LR/PACBIO_DATA/TIPS_212256257_12Mo_LR.hifi_reads.bam \
                prefix=asthma_19 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_19 -j 64 --use-conda" \
                --job-name=asthma_19 
            

            #### number 72
            sbatch --partition standard --wrap "snakemake  -s find_MGE.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_212256257_1Mo_LR/PACBIO_DATA/TIPS_212256257_1Mo_LR.hifi_reads.bam \
                prefix=asthma_20 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_20 -j 64 --use-conda" \
                --job-name=asthma_20 
            
