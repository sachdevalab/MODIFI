

            #### number 2
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208058140_1Mo_LR/PACBIO_DATA/TIPS_208058140_1Mo_LR.hifi_reads.bam \
                prefix=asthma_2 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_2 -j 64" \
                --job-name=asthma_2
            

            #### number 3
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208062141_12Mo_LR/PACBIO_DATA/TIPS_208062141_12Mo_LR.hifi_reads.bam \
                prefix=asthma_3 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_3 -j 64" \
                --job-name=asthma_3
            

            #### number 4
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208062141_1Mo_LR/PACBIO_DATA/TIPS_208062141_1Mo_LR.hifi_reads.bam \
                prefix=asthma_4 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_4 -j 64" \
                --job-name=asthma_4
            

            #### number 5
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208075159_12Mo_LR/PACBIO_DATA/TIPS_208075159_12Mo_LR.hifi_reads.bam \
                prefix=asthma_5 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_5 -j 64" \
                --job-name=asthma_5
            

            #### number 6
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208075159_1Mo_LR/PACBIO_DATA/TIPS_208075159_1Mo_LR.hifi_reads.bam \
                prefix=asthma_6 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_6 -j 64" \
                --job-name=asthma_6
            

            #### number 7
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208107157_12Mo_LR/PACBIO_DATA/TIPS_208107157_12Mo_LR.hifi_reads.bam \
                prefix=asthma_7 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_7 -j 64" \
                --job-name=asthma_7
            

            #### number 8
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_208107157_1Mo_LR/PACBIO_DATA/TIPS_208107157_1Mo_LR.hifi_reads.bam \
                prefix=asthma_8 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_8 -j 64" \
                --job-name=asthma_8
            

            #### number 9
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209129175_12Mo_LR/PACBIO_DATA/TIPS_209129175_12Mo_LR.hifi_reads.bam \
                prefix=asthma_9 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_9 -j 64" \
                --job-name=asthma_9
            

            #### number 10
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209129175_1Mo_LR/PACBIO_DATA/TIPS_209129175_1Mo_LR.hifi_reads.bam \
                prefix=asthma_10 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_10 -j 64" \
                --job-name=asthma_10
            

            #### number 11
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209134177_12Mo_LR/PACBIO_DATA/TIPS_209134177_12Mo_LR.hifi_reads.bam \
                prefix=asthma_11 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_11 -j 64" \
                --job-name=asthma_11
            

            #### number 12
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209134177_1Mo_LR/PACBIO_DATA/TIPS_209134177_1Mo_LR.hifi_reads.bam \
                prefix=asthma_12 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_12 -j 64" \
                --job-name=asthma_12
            

            #### number 13
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209158193_12Mo_LR/PACBIO_DATA/TIPS_209158193_12Mo_LR.hifi_reads.bam \
                prefix=asthma_13 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_13 -j 64" \
                --job-name=asthma_13
            

            #### number 14
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209158193_1Mo_LR/PACBIO_DATA/TIPS_209158193_1Mo_LR.hifi_reads.bam \
                prefix=asthma_14 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_14 -j 64" \
                --job-name=asthma_14
            

            #### number 15
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_210173216_12Mo_LR/PACBIO_DATA/TIPS_210173216_12Mo_LR.hifi_reads.bam \
                prefix=asthma_15 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_15 -j 64" \
                --job-name=asthma_15
            

            #### number 16
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_210173216_1Mo_LR/PACBIO_DATA/TIPS_210173216_1Mo_LR.hifi_reads.bam \
                prefix=asthma_16 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_16 -j 64" \
                --job-name=asthma_16
            

            #### number 17
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_211245247_12Mo_LR/PACBIO_DATA/TIPS_211245247_12Mo_LR.hifi_reads.bam \
                prefix=asthma_17 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_17 -j 64" \
                --job-name=asthma_17
            

            #### number 18
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_211245247_1Mo_LR/PACBIO_DATA/TIPS_211245247_1Mo_LR.hifi_reads.bam \
                prefix=asthma_18 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_18 -j 64" \
                --job-name=asthma_18
            

            #### number 19
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_212256257_12Mo_LR/PACBIO_DATA/TIPS_212256257_12Mo_LR.hifi_reads.bam \
                prefix=asthma_19 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_19 -j 64" \
                --job-name=asthma_19
            

            #### number 20
            sbatch --partition standard --wrap "snakemake -s assembly.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_212256257_1Mo_LR/PACBIO_DATA/TIPS_212256257_1Mo_LR.hifi_reads.bam \
                prefix=asthma_20 \
                work_dir=/groups/banfield/projects/multienv/methylation_temp/asthma_run/asthma_20 -j 64" \
                --job-name=asthma_20
            
