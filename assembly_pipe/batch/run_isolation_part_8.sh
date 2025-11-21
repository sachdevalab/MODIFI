snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR13008124.ccs.bam \
                        prefix=SRR13008124 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR13008124 \
                        -j 64  
