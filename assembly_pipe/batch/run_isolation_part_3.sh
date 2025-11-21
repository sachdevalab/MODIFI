snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR13342926.ccs.bam \
                        prefix=ERR13342926 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR13342926 \
                        -j 64  
