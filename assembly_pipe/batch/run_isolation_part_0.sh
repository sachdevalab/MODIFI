snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR16962737.ccs.bam \
                        prefix=SRR16962737 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR16962737 \
                        -j 64  
