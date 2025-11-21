snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR12778246.ccs.bam \
                        prefix=ERR12778246 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR12778246 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR9682697.ccs.bam \
                        prefix=ERR9682697 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR9682697 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR31303466.ccs.bam \
                        prefix=SRR31303466 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR31303466 \
                        -j 64  
