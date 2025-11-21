snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR10357025.ccs.bam \
                        prefix=ERR10357025 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR10357025 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR6797423.ccs.bam \
                        prefix=ERR6797423 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR6797423 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR16962737.ccs.bam \
                        prefix=SRR16962737 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR16962737 \
                        -j 64  
