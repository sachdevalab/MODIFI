snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR13342922.ccs.bam \
                        prefix=ERR13342922 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR13342922 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/ERR9682698.ccs.bam \
                        prefix=ERR9682698 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//ERR9682698 \
                        -j 64  
snakemake -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR32364866.ccs.bam \
                        prefix=SRR32364866 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR32364866 \
                        -j 64  
