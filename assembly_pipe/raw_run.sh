
            #### number 1
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/RuReacBro_20230708_11_72h_200ppm_r1_LR.hifi_reads.bam \
                prefix=RuReacBro_20230708_11_72h_200ppm_r1_LR \
                work_dir=/home/shuaiw/borg/paper/run/RuReacBro_20230708_11_72h_200ppm_r1_LR" \
                --job-name=RuReacBro_20230708_11_72h_200ppm_r1_LR
            

            #### number 2
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/RuReacBro_20230708_26_72h_NC_r4_LR.hifi_reads.bam \
                prefix=RuReacBro_20230708_26_72h_NC_r4_LR \
                work_dir=/home/shuaiw/borg/paper/run/RuReacBro_20230708_26_72h_NC_r4_LR" \
                --job-name=RuReacBro_20230708_26_72h_NC_r4_LR
            

            #### number 3
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/XSDIA_20240410_R84050_PL8707-001_1-1-A01_bc2007-bc2007.hifi_reads.bam \
                prefix=XSDIA_20240410_R84050_PL8707-001_1-1-A01_bc2007-bc2007 \
                work_dir=/home/shuaiw/borg/paper/run/XSDIA_20240410_R84050_PL8707-001_1-1-A01_bc2007-bc2007" \
                --job-name=XSDIA_20240410_R84050_PL8707-001_1-1-A01_bc2007-bc2007
            

            #### number 4
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/XSDIA_20240410_R84050_PL8708-001_1-1-A01_bc2008-bc2008.hifi_reads.bam \
                prefix=XSDIA_20240410_R84050_PL8708-001_1-1-A01_bc2008-bc2008 \
                work_dir=/home/shuaiw/borg/paper/run/XSDIA_20240410_R84050_PL8708-001_1-1-A01_bc2008-bc2008" \
                --job-name=XSDIA_20240410_R84050_PL8708-001_1-1-A01_bc2008-bc2008
            

            #### number 5
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow_bioreactor/XSDIA_20240410_R84050_PL8709-001_1-1-A01_bc2096-bc2096.hifi_reads.bam \
                prefix=XSDIA_20240410_R84050_PL8709-001_1-1-A01_bc2096-bc2096 \
                work_dir=/home/shuaiw/borg/paper/run/XSDIA_20240410_R84050_PL8709-001_1-1-A01_bc2096-bc2096" \
                --job-name=XSDIA_20240410_R84050_PL8709-001_1-1-A01_bc2096-bc2096
            

            #### number 6
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/cow/RuReacBro_20230708_Comb_RF_LR.hifi_reads.bam \
                prefix=RuReacBro_20230708_Comb_RF_LR \
                work_dir=/home/shuaiw/borg/paper/run/RuReacBro_20230708_Comb_RF_LR" \
                --job-name=RuReacBro_20230708_Comb_RF_LR
            

            #### number 7
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR12723528/ERR12723528.ccs.bam \
                prefix=ERR12723528 \
                work_dir=/home/shuaiw/borg/paper/run/ERR12723528" \
                --job-name=ERR12723528
            

            #### number 8
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR12723529/ERR12723529.ccs.bam \
                prefix=ERR12723529 \
                work_dir=/home/shuaiw/borg/paper/run/ERR12723529" \
                --job-name=ERR12723529
            

            #### number 9   too less data
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621427/ERR5621427.ccs.bam \
                prefix=ERR5621427 \
                work_dir=/home/shuaiw/borg/paper/run/ERR5621427" \
                --job-name=ERR5621427
            

            #### number 10
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621429/ERR5621429.ccs.bam \
                prefix=ERR5621429 \
                work_dir=/home/shuaiw/borg/paper/run/ERR5621429" \
                --job-name=ERR5621429
            

            #### number 11
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ERR5621430/ERR5621430.ccs.bam \
                prefix=ERR5621430 \
                work_dir=/home/shuaiw/borg/paper/run/ERR5621430" \
                --job-name=ERR5621430
            

            #### number 12
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240037_3G1_pacbio.bam \
                prefix=NANO_2_INF1240037_3G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1240037_3G1_pacbio" \
                --job-name=NANO_2_INF1240037_3G1_pacbio
            

            #### number 13
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_2G1_pacbio.bam \
                prefix=NANO_2_INF1240040_2G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1240040_2G1_pacbio" \
                --job-name=NANO_2_INF1240040_2G1_pacbio
            

            #### number 14
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_3G1_pacbio.bam \
                prefix=NANO_2_INF1240040_3G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1240040_3G1_pacbio" \
                --job-name=NANO_2_INF1240040_3G1_pacbio
            

            #### number 15
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_4G1_pacbio.bam \
                prefix=NANO_2_INF1240040_4G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1240040_4G1_pacbio" \
                --job-name=NANO_2_INF1240040_4G1_pacbio
            

            #### number 16
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1310010_2G1_pacbio.bam \
                prefix=NANO_2_INF1310010_2G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1310010_2G1_pacbio" \
                --job-name=NANO_2_INF1310010_2G1_pacbio
            

            #### number 17
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1310010_3G1_pacbio.bam \
                prefix=NANO_2_INF1310010_3G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1310010_3G1_pacbio" \
                --job-name=NANO_2_INF1310010_3G1_pacbio
            

            #### number 18
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330001_4G1_pacbio.bam \
                prefix=NANO_2_INF1330001_4G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1330001_4G1_pacbio" \
                --job-name=NANO_2_INF1330001_4G1_pacbio
            

            #### number 19
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330004_3G1_pacbio.bam \
                prefix=NANO_2_INF1330004_3G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1330004_3G1_pacbio" \
                --job-name=NANO_2_INF1330004_3G1_pacbio
            

            #### number 20
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1330004_4G1_pacbio.bam \
                prefix=NANO_2_INF1330004_4G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1330004_4G1_pacbio" \
                --job-name=NANO_2_INF1330004_4G1_pacbio
            

            #### number 21
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_2G1_pacbio.bam \
                prefix=NANO_2_INF1340011_2G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1340011_2G1_pacbio" \
                --job-name=NANO_2_INF1340011_2G1_pacbio
            

            #### number 22
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_3G1_pacbio.bam \
                prefix=NANO_2_INF1340011_3G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1340011_3G1_pacbio" \
                --job-name=NANO_2_INF1340011_3G1_pacbio
            

            #### number 23
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1340011_4G1_pacbio.bam \
                prefix=NANO_2_INF1340011_4G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_INF1340011_4G1_pacbio" \
                --job-name=NANO_2_INF1340011_4G1_pacbio
            

            #### number 24
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_MAT1240037_2G1_pacbio.bam \
                prefix=NANO_2_MAT1240037_2G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_2_MAT1240037_2G1_pacbio" \
                --job-name=NANO_2_MAT1240037_2G1_pacbio
            

            #### number 25
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_5G1_pacbio.bam \
                prefix=NANO_3_INF1240040_5G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1240040_5G1_pacbio" \
                --job-name=NANO_3_INF1240040_5G1_pacbio
            

            #### number 26
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1240040_6G1_pacbio.bam \
                prefix=NANO_3_INF1240040_6G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1240040_6G1_pacbio" \
                --job-name=NANO_3_INF1240040_6G1_pacbio
            

            #### number 27
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1310001_8G1_pacbio.bam \
                prefix=NANO_3_INF1310001_8G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1310001_8G1_pacbio" \
                --job-name=NANO_3_INF1310001_8G1_pacbio
            

            #### number 28
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1310007_7G1_pacbio.bam \
                prefix=NANO_3_INF1310007_7G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1310007_7G1_pacbio" \
                --job-name=NANO_3_INF1310007_7G1_pacbio
            

            #### number 29
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_5G1_pacbio.bam \
                prefix=NANO_3_INF1330004_5G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1330004_5G1_pacbio" \
                --job-name=NANO_3_INF1330004_5G1_pacbio
            

            #### number 30
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_7G1_pacbio.bam \
                prefix=NANO_3_INF1330004_7G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1330004_7G1_pacbio" \
                --job-name=NANO_3_INF1330004_7G1_pacbio
            

            #### number 31
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1330004_8G1_pacbio.bam \
                prefix=NANO_3_INF1330004_8G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1330004_8G1_pacbio" \
                --job-name=NANO_3_INF1330004_8G1_pacbio
            

            #### number 32
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340008_8G1_pacbio.bam \
                prefix=NANO_3_INF1340008_8G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1340008_8G1_pacbio" \
                --job-name=NANO_3_INF1340008_8G1_pacbio
            

            #### number 33
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340011_6G1_pacbio.bam \
                prefix=NANO_3_INF1340011_6G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1340011_6G1_pacbio" \
                --job-name=NANO_3_INF1340011_6G1_pacbio
            

            #### number 34
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_3_INF1340011_7G1_pacbio.bam \
                prefix=NANO_3_INF1340011_7G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_3_INF1340011_7G1_pacbio" \
                --job-name=NANO_3_INF1340011_7G1_pacbio
            

            #### number 35
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240043_5G1_pacbio.bam \
                prefix=NANO_4_INF1240043_5G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_4_INF1240043_5G1_pacbio" \
                --job-name=NANO_4_INF1240043_5G1_pacbio
            

            #### number 36
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240045_6G1_pacbio.bam \
                prefix=NANO_4_INF1240045_6G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_4_INF1240045_6G1_pacbio" \
                --job-name=NANO_4_INF1240045_6G1_pacbio
            

            #### number 37
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1240119_7G1_pacbio.bam \
                prefix=NANO_4_INF1240119_7G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_4_INF1240119_7G1_pacbio" \
                --job-name=NANO_4_INF1240119_7G1_pacbio
            

            #### number 38
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1330101_8G1_pacbio.bam \
                prefix=NANO_4_INF1330101_8G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_4_INF1330101_8G1_pacbio" \
                --job-name=NANO_4_INF1330101_8G1_pacbio
            

            #### number 39
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/infant/NANO_4_INF1340021_5G1_pacbio.bam \
                prefix=NANO_4_INF1340021_5G1_pacbio \
                work_dir=/home/shuaiw/borg/paper/run/NANO_4_INF1340021_5G1_pacbio" \
                --job-name=NANO_4_INF1340021_5G1_pacbio
            

            #### number 40
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/ocean/pbio-2857.29455.bc1003_BAK8A_OA--bc1003_BAK8A_OA.hifi_reads.bc1003_BAK8A_OA.ccs.bam \
                prefix=pbio-2857 \
                work_dir=/home/shuaiw/borg/paper/run/pbio-2857 --rerun-incomplete" \
                --job-name=pbio-2857
            

            #### number 41
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2021.bam \
                prefix=m84039_230624_013044_s3 \
                work_dir=/home/shuaiw/borg/paper/run/m84039_230624_013044_s3" \
                --job-name=m84039_230624_013044_s3
            

            #### number 42
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230624_013044_s3.hifi_reads.bc2023.bam \
                prefix=m84039_230624_013044_s3 \
                work_dir=/home/shuaiw/borg/paper/run/m84039_230624_013044_s3" \
                --job-name=m84039_230624_013044_s3
            

            #### number 43
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2022.bam \
                prefix=m84039_230626_214113_s4 \
                work_dir=/home/shuaiw/borg/paper/run/m84039_230626_214113_s4" \
                --job-name=m84039_230626_214113_s4
            

            #### number 44
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_214113_s4.hifi_reads.bc2024.bam \
                prefix=m84039_230626_214113_s4 \
                work_dir=/home/shuaiw/borg/paper/run/m84039_230626_214113_s4" \
                --job-name=m84039_230626_214113_s4
            

            #### number 45
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2025.bam \
                prefix=m84039_230626_221130_s1 \
                work_dir=/home/shuaiw/borg/paper/run/m84039_230626_221130_s1" \
                --job-name=m84039_230626_221130_s1
            

            #### number 46
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/m84039_230626_221130_s1.hifi_reads.bc2026.bam \
                prefix=m84039_230626_221130_s1 \
                work_dir=/home/shuaiw/borg/paper/run/m84039_230626_221130_s1" \
                --job-name=m84039_230626_221130_s1
            

            #### number 47
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/soil/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.bam \
                prefix=XRSBK_20221007_S64018_PL100268287-1_C01 \
                work_dir=/home/shuaiw/borg/paper/run/XRSBK_20221007_S64018_PL100268287-1_C01" \
                --job-name=XRSBK_20221007_S64018_PL100268287-1_C01
            

            #### number 48
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR14074352/SRR14074352.ccs.bam \
                prefix=SRR14074352 \
                work_dir=/home/shuaiw/borg/paper/run/SRR14074352" \
                --job-name=SRR14074352
            

            #### number 49
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR23446539/SRR23446539.ccs.bam \
                prefix=SRR23446539 \
                work_dir=/home/shuaiw/borg/paper/run/SRR23446539" \
                --job-name=SRR23446539
            

            #### number 50
            sbatch --partition standard --wrap "snakemake --config \
                hifi_bam=/home/shuaiw/borg/paper/aws/SRR23446540/SRR23446540.ccs.bam \
                prefix=SRR23446540 \
                work_dir=/home/shuaiw/borg/paper/run/SRR23446540" \
                --job-name=SRR23446540


sbatch --partition standard --wrap "snakemake -s methylation.smk --config                 hifi_bam=/home/shuaiw/borg/paper/aws/ocean/pbio-2857.29455.bc1003_BAK8A_OA--bc1003_BAK8A_OA.hifi_reads.bc1003_BAK8A_OA.ccs.bam                 prefix=pbio-2857                 work_dir=/home/shuaiw/borg/paper/run/pbio-2857 --rerun-incomplete" --job-name=pbio-2857
            
