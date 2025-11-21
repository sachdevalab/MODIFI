snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR6536208.ccs.bam \
                        prefix=ERR6536208 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR6536208 \
                        -j 64  
