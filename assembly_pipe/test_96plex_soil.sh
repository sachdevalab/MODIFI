sbatch --partition standard --wrap "snakemake --config \
    hifi_bam=/home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.raw.p05.p05.bam \
    prefix=meta_p05 \
    work_dir=/home/shuaiw/borg/paper/test/meta_p05 -j 64" \
    --job-name=meta_p05
