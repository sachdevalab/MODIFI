snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR27457941.ccs.bam \
                        prefix=SRR27457941 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR27457941 \
                        -j 64  
