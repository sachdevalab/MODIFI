snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR31014709.ccs.bam \
                        prefix=SRR31014709 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR31014709 \
                        -j 64  
