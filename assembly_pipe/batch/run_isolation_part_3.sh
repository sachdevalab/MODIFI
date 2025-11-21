snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR13342925.ccs.bam \
                        prefix=ERR13342925 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR13342925 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR9682725.ccs.bam \
                        prefix=ERR9682725 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR9682725 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR32364929.ccs.bam \
                        prefix=SRR32364929 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR32364929 \
                        -j 64  
