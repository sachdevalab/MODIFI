sbatch --partition standard --wrap "snakemake --config \
    hifi_bam=/home/shuaiw/borg/paper/linkage/meta/m64004_210929_143746.raw.p50.soil.merge.bam \
    prefix=meta_p50 \
    work_dir=/home/shuaiw/borg/paper/test/meta_p50 -j 64" \
    --job-name=meta_p50
