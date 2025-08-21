
            #### number 1
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.hifi_reads.bam \
                prefix=96plex \
                work_dir=/home/shuaiw/borg/paper/run2/96plex -j 64" \
                --job-name=96plex
            

            #### number 2
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/RuReacBro_20230708_11_72h_200ppm_r1_LR.hifi_reads.bam \
                prefix=cow_bioreactor_1 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_bioreactor_1 -j 64" \
                --job-name=cow_bioreactor_1
            

            #### number 3
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/RuReacBro_20230708_26_72h_NC_r4_LR.hifi_reads.bam \
                prefix=cow_bioreactor_2 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_bioreactor_2 -j 64" \
                --job-name=cow_bioreactor_2
            

            #### number 4
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/XSDIA_20240410_R84050_PL8707-001_1-1-A01_bc2007-bc2007.hifi_reads.bam \
                prefix=cow_bioreactor_3 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_bioreactor_3 -j 64" \
                --job-name=cow_bioreactor_3
            

            #### number 5
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/XSDIA_20240410_R84050_PL8708-001_1-1-A01_bc2008-bc2008.hifi_reads.bam \
                prefix=cow_bioreactor_4 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_bioreactor_4 -j 64" \
                --job-name=cow_bioreactor_4
            

            #### number 6
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/XSDIA_20240410_R84050_PL8709-001_1-1-A01_bc2096-bc2096.hifi_reads.bam \
                prefix=cow_bioreactor_5 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_bioreactor_5 -j 64" \
                --job-name=cow_bioreactor_5
            

            #### number 7
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow/RuReacBro_20230708_Comb_RF_LR.hifi_reads.bam \
                prefix=cow_1 \
                work_dir=/home/shuaiw/borg/paper/run2/cow_1 -j 64" \
                --job-name=cow_1
            

            #### number 8
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR12723528/ERR12723528.ccs.bam \
                prefix=ERR12723528_mice \
                work_dir=/home/shuaiw/borg/paper/run2/ERR12723528_mice -j 64" \
                --job-name=ERR12723528_mice
            

            #### number 9
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR12723529/ERR12723529.ccs.bam \
                prefix=ERR12723529_mice \
                work_dir=/home/shuaiw/borg/paper/run2/ERR12723529_mice -j 64" \
                --job-name=ERR12723529_mice
            

            #### number 10
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621427/ERR5621427.ccs.bam \
                prefix=ERR5621427_sludge \
                work_dir=/home/shuaiw/borg/paper/run2/ERR5621427_sludge -j 64" \
                --job-name=ERR5621427_sludge
            

            #### number 11
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621429/ERR5621429.ccs.bam \
                prefix=ERR5621429_sludge \
                work_dir=/home/shuaiw/borg/paper/run2/ERR5621429_sludge -j 64" \
                --job-name=ERR5621429_sludge
            

            #### number 12
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621430/ERR5621430.ccs.bam \
                prefix=ERR5621430_sludge \
                work_dir=/home/shuaiw/borg/paper/run2/ERR5621430_sludge -j 64" \
                --job-name=ERR5621430_sludge
            

            #### number 13
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240037_3G1_pacbio.bam \
                prefix=infant_1 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_1 -j 64" \
                --job-name=infant_1
            

            #### number 14
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_2G1_pacbio.bam \
                prefix=infant_2 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_2 -j 64" \
                --job-name=infant_2
            

            #### number 15
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_3G1_pacbio.bam \
                prefix=infant_3 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_3 -j 64" \
                --job-name=infant_3
            

            #### number 16
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_4G1_pacbio.bam \
                prefix=infant_4 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_4 -j 64" \
                --job-name=infant_4
            

            #### number 17
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1310010_2G1_pacbio.bam \
                prefix=infant_5 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_5 -j 64" \
                --job-name=infant_5
            

            #### number 18
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1310010_3G1_pacbio.bam \
                prefix=infant_6 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_6 -j 64" \
                --job-name=infant_6
            

            #### number 19
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330001_4G1_pacbio.bam \
                prefix=infant_7 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_7 -j 64" \
                --job-name=infant_7
            

            #### number 20
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330004_3G1_pacbio.bam \
                prefix=infant_8 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_8 -j 64" \
                --job-name=infant_8
            

            #### number 21
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330004_4G1_pacbio.bam \
                prefix=infant_9 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_9 -j 64" \
                --job-name=infant_9
            

            #### number 22
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_2G1_pacbio.bam \
                prefix=infant_10 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_10 -j 64" \
                --job-name=infant_10
            

            #### number 23
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_3G1_pacbio.bam \
                prefix=infant_11 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_11 -j 64" \
                --job-name=infant_11
            

            #### number 24
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_4G1_pacbio.bam \
                prefix=infant_12 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_12 -j 64" \
                --job-name=infant_12
            

            #### number 25
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_MAT1240037_2G1_pacbio.bam \
                prefix=infant_13 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_13 -j 64" \
                --job-name=infant_13
            

            #### number 26
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_5G1_pacbio.bam \
                prefix=infant_14 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_14 -j 64" \
                --job-name=infant_14
            

            #### number 27
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_6G1_pacbio.bam \
                prefix=infant_15 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_15 -j 64" \
                --job-name=infant_15
            

            #### number 28
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1310001_8G1_pacbio.bam \
                prefix=infant_16 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_16 -j 64" \
                --job-name=infant_16
            

            #### number 29
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1310007_7G1_pacbio.bam \
                prefix=infant_17 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_17 -j 64" \
                --job-name=infant_17
            

            #### number 30
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_5G1_pacbio.bam \
                prefix=infant_18 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_18 -j 64" \
                --job-name=infant_18
            

            #### number 31
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_7G1_pacbio.bam \
                prefix=infant_19 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_19 -j 64" \
                --job-name=infant_19
            

            #### number 32
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_8G1_pacbio.bam \
                prefix=infant_20 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_20 -j 64" \
                --job-name=infant_20
            

            #### number 33
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340008_8G1_pacbio.bam \
                prefix=infant_21 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_21 -j 64" \
                --job-name=infant_21
            

            #### number 34
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340011_6G1_pacbio.bam \
                prefix=infant_22 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_22 -j 64" \
                --job-name=infant_22
            

            #### number 35
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340011_7G1_pacbio.bam \
                prefix=infant_23 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_23 -j 64" \
                --job-name=infant_23
            

            #### number 36
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240043_5G1_pacbio.bam \
                prefix=infant_24 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_24 -j 64" \
                --job-name=infant_24
            

            #### number 37
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240045_6G1_pacbio.bam \
                prefix=infant_25 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_25 -j 64" \
                --job-name=infant_25
            

            #### number 38
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240119_7G1_pacbio.bam \
                prefix=infant_26 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_26 -j 64" \
                --job-name=infant_26
            

            #### number 39
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1330101_8G1_pacbio.bam \
                prefix=infant_27 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_27 -j 64" \
                --job-name=infant_27
            

            #### number 40
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1340021_5G1_pacbio.bam \
                prefix=infant_28 \
                work_dir=/home/shuaiw/borg/paper/run2/infant_28 -j 64" \
                --job-name=infant_28
            

            #### number 41
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ocean/pbio-2857.29455.bc1003_BAK8A_OA--bc1003_BAK8A_OA.hifi_reads.bc1003_BAK8A_OA.ccs.bam \
                prefix=ocean_1 \
                work_dir=/home/shuaiw/borg/paper/run2/ocean_1 -j 64" \
                --job-name=ocean_1
            

            #### number 42
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
                prefix=soil_s3_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s3_1 -j 64" \
                --job-name=soil_s3_1
            

            #### number 43
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                prefix=soil_s3_2 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s3_2 -j 64" \
                --job-name=soil_s3_2
            

            #### number 44
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
                prefix=soil_s4_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s4_1 -j 64" \
                --job-name=soil_s4_1
            

            #### number 45
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2024.bam \
                prefix=soil_s4_2 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s4_2 -j 64" \
                --job-name=soil_s4_2
            

            #### number 46
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2025.bam \
                prefix=soil_s1_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s1_1 -j 64" \
                --job-name=soil_s1_1
            

            #### number 47
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
                prefix=soil_s1_2 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_s1_2 -j 64" \
                --job-name=soil_s1_2
            

            #### number 48
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam \
                prefix=soil_1 \
                work_dir=/home/shuaiw/borg/paper/run2/soil_1 -j 64" \
                --job-name=soil_1
            

            #### number 49
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR14074352/SRR14074352.ccs.bam \
                prefix=SRR14074352_human \
                work_dir=/home/shuaiw/borg/paper/run2/SRR14074352_human -j 64" \
                --job-name=SRR14074352_human
            

            #### number 50
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR23446539/SRR23446539.ccs.bam \
                prefix=SRR23446539_sugarcane \
                work_dir=/home/shuaiw/borg/paper/run2/SRR23446539_sugarcane -j 64" \
                --job-name=SRR23446539_sugarcane
            

            #### number 51
            sbatch --partition standard --wrap "snakemake -s drep.smk --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR23446540/SRR23446540.ccs.bam \
                prefix=SRR23446540_sugarcane \
                work_dir=/home/shuaiw/borg/paper/run2/SRR23446540_sugarcane -j 64" \
                --job-name=SRR23446540_sugarcane
            
