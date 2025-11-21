snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR13342926.ccs.bam \
                        prefix=ERR13342926 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR13342926 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR9682728.ccs.bam \
                        prefix=ERR9682728 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR9682728 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR32375313.ccs.bam \
                        prefix=SRR32375313 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR32375313 \
                        -j 64  
