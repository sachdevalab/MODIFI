snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR6797423.ccs.bam \
                        prefix=ERR6797423 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR6797423 \
                        -j 64  
