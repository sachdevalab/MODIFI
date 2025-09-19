# snakemake -s isolation.smk \
#     --config hifi_bam=/home/shuaiw/borg/bench/soil_zymo/run4/E_coli_H10407.bam \
#     prefix=test2 \
#     work_dir=/home/shuaiw/borg/paper/isolation/test2 \
#     -j 64


 sbatch  --partition standard --wrap "snakemake -s isolation.smk \
    --config hifi_bam=/home/shuaiw/borg/paper/aws/isolate/archaea/SRR31014709/SRR31014709.ccs.bam \
    prefix=SRR31014709 \
    work_dir=/home/shuaiw/borg/paper/isolation/archaea/SRR31014709 \
    -j 64"  --job-name=iso_1

 sbatch  --partition standard --wrap "snakemake -s isolation.smk \
    --config hifi_bam=/home/shuaiw/borg/paper/aws/isolate/archaea/SRR27457941/SRR27457941.ccs.bam \
    prefix=SRR27457941 \
    work_dir=/home/shuaiw/borg/paper/isolation/archaea/SRR27457941 \
    -j 64"  --job-name=iso_2

