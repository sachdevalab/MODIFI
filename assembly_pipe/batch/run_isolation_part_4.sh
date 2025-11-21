snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR13655791.ccs.bam \
                        prefix=ERR13655791 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR13655791 \
                        -j 64  
