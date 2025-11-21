snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR6536208.ccs.bam \
                        prefix=ERR6536208 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR6536208 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR15032816.ccs.bam \
                        prefix=SRR15032816 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR15032816 \
                        -j 64  
