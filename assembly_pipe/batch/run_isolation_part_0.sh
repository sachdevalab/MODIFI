snakemake --rerun-incomplete -s isolation.smk \
                        --config hifi_bam=/home/shuaiw/borg/paper/isolation/batch2_ccs_bam/SRR19519931.ccs.bam \
                        prefix=SRR19519931 \
                        work_dir=/home/shuaiw/borg/paper/isolation/batch2_results//SRR19519931 \
                        -j 64  
