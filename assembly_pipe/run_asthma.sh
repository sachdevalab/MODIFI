
            # #### number 9
            snakemake -s annotation.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209129175_12Mo_LR/PACBIO_DATA/TIPS_209129175_12Mo_LR.hifi_reads.bam \
                prefix=asthma_9 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_9 -j 64 
            

            # #### number 14
            snakemake -s annotation.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_209158193_1Mo_LR/PACBIO_DATA/TIPS_209158193_1Mo_LR.hifi_reads.bam \
                prefix=asthma_14 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_14 -j 64 
            

            # #### number 18
            snakemake -s annotation.smk --config \
                hifi_bam=/groups/diamond/sequences/2024/TIPS_211245247_1Mo_LR/PACBIO_DATA/TIPS_211245247_1Mo_LR.hifi_reads.bam \
                prefix=asthma_18 \
                work_dir=/home/shuaiw/borg/paper/run2/asthma_18 -j 64 
            
