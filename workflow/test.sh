sbatch --partition standard --wrap "/usr/bin/time -v -o large.time snakemake --config whole_bam=/home/shuaiw/methylation/data/borg/large_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam whole_ref=/home/shuaiw/methylation/data/borg/contigs/large.contigs.fa work_dir=/home/shuaiw/methylation/data/borg/new_test6"
sbatch --partition standard --wrap "/usr/bin/time -v -o test_100.time snakemake --config whole_bam=/home/shuaiw/borg/test_100/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam whole_ref=/home/shuaiw/borg/contigs/test_100.fa work_dir=/home/shuaiw/methylation/data/borg/new_test7"
snakemake --config whole_bam=/home/shuaiw/methylation/data/borg/new_test7/bams/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_91_L.bam whole_ref=/home/shuaiw/methylation/data/borg/new_test7/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_91_L.fa work_dir=/home/shuaiw/methylation/data/borg/new_test8

sbatch --partition standard --wrap "/usr/bin/time -v -o large.time snakemake --config whole_bam=/home/shuaiw/borg/seven_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam  whole_ref=/home/shuaiw/borg/contigs/seven_contigs.fasta  work_dir=/home/shuaiw/methylation/data/borg/new_test10" --job-name=7
sbatch --partition standard --wrap "/usr/bin/time -v -o test_200.time snakemake --config whole_bam=/home/shuaiw/methylation/data/borg/test_200/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam whole_ref=/home/shuaiw/borg/contigs/test_200.fa work_dir=/home/shuaiw/methylation/data/borg/test_200_our" --job-name=test_200 


sbatch --partition standard --wrap "/usr/bin/time -v -o test_500.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/test_500/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
whole_ref=/home/shuaiw/borg/contigs/test_500.fa\
 work_dir=/home/shuaiw/methylation/data/borg/test_500_our_3"\
  --job-name=test_500

sbatch --partition standard --wrap "/usr/bin/time -v -o test_300.time snakemake --config whole_bam=/
home/shuaiw/methylation/data/borg/test_300/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam whole_ref=/home/shuaiw/borg/contigs/test_300.fa w
ork_dir=/home/shuaiw/methylation/data/borg/test_300_our" --job-name=test_300


sbatch --partition standard --wrap "/usr/bin/time -v -o all.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa\
 work_dir=/home/shuaiw/methylation/data/borg/all_test"\
  --job-name=all


sbatch --partition standard --wrap "/usr/bin/time -v -o seven.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/seven_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/seven_contigs.fasta\
 work_dir=/home/shuaiw/methylation/data/borg/seven_test"\
  --job-name=seven

sbatch --partition standard --wrap "/usr/bin/time -v -o SAMN07447446.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/rebase/SAMN07447446/align2.bam \
whole_ref=/home/shuaiw/methylation/data/rebase/SAMN07447446/meta.fa\
 work_dir=/home/shuaiw/methylation/data/rebase/SAMN07447446/result2"\
  --job-name=SAMN07447446



sbatch --partition standard --wrap "/usr/bin/time -v -o human.time snakemake \
--config whole_bam=/home/shuaiw/borg/human/human_000733.subreads.align.bam \
whole_ref=/home/shuaiw/borg/hg38/GCF_000001405.40_GRCh38.p14_genomic.fasta\
 work_dir=/home/shuaiw/borg/human_test"\
  --job-name=human

  snakemake --config whole_bam=/home/shuaiw/methylation/data/borg/b_contigs/11.align.bam\
 whole_ref=/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa\
  work_dir=/home/shuaiw/methylation/data/borg/new_test10


  snakemake -s pipeline2.smk --config whole_bam=/home/shuaiw/methylation/data/borg/b_contigs/11.align.bam \
  whole_ref=/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa \
  work_dir=/home/shuaiw/methylation/data/borg/new_test11