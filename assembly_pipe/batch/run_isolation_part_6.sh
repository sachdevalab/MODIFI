snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR13655871.ccs.bam \
                        prefix=ERR13655871 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR13655871 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR13008124.ccs.bam \
                        prefix=SRR13008124 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR13008124 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR34878266.ccs.bam \
                        prefix=SRR34878266 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR34878266 \
                        -j 64  
