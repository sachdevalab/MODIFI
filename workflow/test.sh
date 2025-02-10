sbatch --partition standard --wrap "/usr/bin/time -v -o large.time snakemake --config whole_bam=/home/shuaiw/methylation/data/borg/large_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam whole_ref=/home/shuaiw/methylation/data/borg/contigs/large.contigs.fa work_dir=/home/shuaiw/methylation/data/borg/new_test6"
sbatch --partition standard --wrap "/usr/bin/time -v -o test_100.time snakemake --config whole_bam=/home/shuaiw/borg/test_100/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam whole_ref=/home/shuaiw/borg/contigs/test_100.fa work_dir=/home/shuaiw/methylation/data/borg/new_test7"

