snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR13655791.ccs.bam \
                        prefix=ERR13655791 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR13655791 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR9682737.ccs.bam \
                        prefix=ERR9682737 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR9682737 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR32918029.ccs.bam \
                        prefix=SRR32918029 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR32918029 \
                        -j 64  
