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
ork_dir=/home/shuaiw/methylation/data/borg/bench/test_300" --job-name=test_300


sbatch --partition standard --wrap "/usr/bin/time -v -o all.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa\
 work_dir=/home/shuaiw/methylation/data/borg/all_test2"\
  --job-name=all


sbatch --partition standard --wrap "/usr/bin/time -v -o seven.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/seven_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/seven_contigs.fasta\
 work_dir=/home/shuaiw/methylation/data/borg/bench/seven_test"\
  --job-name=seven

sbatch --partition standard --wrap "/usr/bin/time -v -o SAMN07447446.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/rebase/SAMN07447446/align2.bam \
whole_ref=/home/shuaiw/methylation/data/rebase/SAMN07447446/meta.fa\
 work_dir=/home/shuaiw/methylation/data/rebase/SAMN07447446/result2"\
  --job-name=SAMN07447446



sbatch --partition standard --wrap "/usr/bin/time -v -o human.time snakemake \
--config whole_bam=/home/shuaiw/borg/human/human_000733.subreads.align.bam \
whole_ref=/home/shuaiw/borg/hg38/GCF_000001405.40_GRCh38.p14_genomic.fasta\
 work_dir=/home/shuaiw/borg/bench/human_test"\
  --job-name=human



   sbatch --partition standard --wrap "/usr/bin/time -v -o borg.time snakemake --config whole_bam=/home/shuaiw/methylation/data/borg/all_borg/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam\
 whole_ref=/home/shuaiw/methylation/data/borg/all_borg/all_borg.fasta\
  work_dir=/home/shuaiw/methylation/data/borg/bench/ccs_borg2_NM0 \
  read_type=ccs max_NM=0" --job-name=borg


  sbatch --partition standard --wrap "/usr/bin/time -v -o all_ccs.time snakemake --rerun-incomplete \
--config whole_bam=/home/shuaiw/methylation/data/borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/all_test_ccs read_type=ccs "\
  --job-name=all_ccs


  sbatch --partition standard --wrap " genomad end-to-end --relaxed --cleanup --enable-score-calibration \
  --threads 64 --sensitivity 7.0 --force-auto \
  /home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
  /home/shuaiw/methylation/data/borg/bench/genomad/ \
  /groups/diamond/databases/genomad/v1.7/" --job-name=genomad


   sbatch --partition standard --wrap "MicrobeMod annotate_rm \
   -f /home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
    -o /home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_RM -t 64" \
   --job-name=MicrobeMod



    sbatch --partition standard --wrap "/usr/bin/time -v -o all_break.time snakemake \
--config whole_bam=/home/shuaiw/borg/all_break/all_break.align.ccs.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/all_break.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/all_break read_type=ccs "\
  --job-name=all_break

  sbatch --partition standard --wrap "/usr/bin/time -v -o all_subreads.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa\
 work_dir=/home/shuaiw/methylation/data/borg/bench/all_subreads"\
  --job-name=all

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_ccs2.time snakemake --rerun-incomplete \
--config whole_bam=/home/shuaiw/methylation/data/borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/all_test_ccs2 read_type=ccs min_len=5000 max_NM=10 min_cov=10"\
  --job-name=all_ccs

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_ccs3.time snakemake --rerun-incomplete \
--config whole_bam=/home/shuaiw/methylation/data/borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/all_test_ccs3 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=all_ccs3



    sbatch --partition standard --wrap "/usr/bin/time -v -o all_s3.time snakemake --rerun-incomplete \
--config whole_bam=/home/shuaiw/borg/all_bams/m84039_230624_013044_s3.hifi_reads.bc2021.align.ccs.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/s3 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=s3

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_s4.time snakemake --rerun-incomplete \
--config whole_bam=/home/shuaiw/borg/all_bams/m84039_230626_214113_s4.hifi_reads.bc2022.align.ccs.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/s4 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=s4

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_s5.time snakemake \
--config whole_bam=/home/shuaiw/borg/all_bams/m84039_230624_013044_s3.hifi_reads.bc2023.align.ccs.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/luis/m84039_230624_013044_s3.hifi_reads.bc2023 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=s5

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_s6.time snakemake \
--config whole_bam=/home/shuaiw/borg/all_bams/m84039_230626_214113_s4.hifi_reads.bc2024.align.ccs.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/luis/m84039_230626_214113_s4.hifi_reads.bc2024 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=s6

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_s7.time snakemake \
--config whole_bam=/home/shuaiw/borg/all_bams/m84039_230626_221130_s1.hifi_reads.bc2025.align.ccs.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/luis/m84039_230626_221130_s1.hifi_reads.bc2025 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=s7

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_s8.time snakemake \
--config whole_bam=/home/shuaiw/borg/all_bams/m84039_230626_221130_s1.hifi_reads.bc2026.align.ccs.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/luis/m84039_230626_221130_s1.hifi_reads.bc2026 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=s8


    sbatch --partition standard --wrap "snakemake \
--config whole_bam=/home/shuaiw/borg/mis_assembly/align.ccs.bam \
whole_ref=/home/shuaiw/borg/contigs/mis_contigs.fasta \
 work_dir=/home/shuaiw/borg/mis_assembly read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=mis

    sbatch --partition standard --wrap "snakemake \
--config whole_bam=/home/shuaiw/borg/allison/NANO_2_INF1330004_4G1.align.ccs.bam \
whole_ref=/home/shuaiw/borg/allison/NANO_2_INF1330004_4PB_HR_HIFIASM_META_scaffold_min1000.fa \
 work_dir=/home/shuaiw/borg/allison/NANO_2_INF1330004_4G1 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=a1

      sbatch --partition standard --wrap "snakemake \
--config whole_bam=/home/shuaiw/borg/allison/NANO_2_INF1340011_4G1.align.ccs.bam \
whole_ref=/home/shuaiw/borg/allison/NANO_2_INF1340011_4PB_HR_HIFIASM_META_scaffold_min1000.fa \
 work_dir=/home/shuaiw/borg/allison/NANO_2_INF1340011_4G1 read_type=ccs min_len=5000 max_NM=30 min_cov=5"\
  --job-name=a2


    sbatch --partition standard --wrap "snakemake \
--config whole_bam=/home/shuaiw/borg/allison/NANO_2_INF1330004_4G1.align.ccs.bam \
whole_ref=/home/shuaiw/borg/allison/NANO_2_INF1330004_4PB_HR_HIFIASM_META_scaffold_min1000.fa \
 work_dir=/home/shuaiw/borg/allison/NANO_2_INF1330004_4G1_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
  --job-name=a1

      sbatch --partition standard --wrap "snakemake \
--config whole_bam=/home/shuaiw/borg/allison/NANO_2_INF1340011_4G1.align.ccs.bam \
whole_ref=/home/shuaiw/borg/allison/NANO_2_INF1340011_4PB_HR_HIFIASM_META_scaffold_min1000.fa \
 work_dir=/home/shuaiw/borg/allison/NANO_2_INF1340011_4G1_NM3 read_type=ccs min_len=5000 max_NM=3 min_cov=5"\
  --job-name=a2

python cal_invasion_score.py --work_dir /home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_more --plasmid_file /groups/diamond/projects/animal/rumen/RuReacBro20203/assembly/RuReacBro_20230708_11_72h_200ppm_r1_LR/genomad/circular_elements_summary/circular_elements_plasmid_summary.tsv

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_ccs_1k.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/all_ccs_1k read_type=ccs min_len=1000 max_NM=10 min_cov=5"\
  --job-name=all_ccs_1k

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_ccs_1k_NM1000.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/all_ccs_1k_NM1000 read_type=ccs min_len=1000 max_NM=1000 min_cov=5 clean=False"\
  --job-name=all_ccs_1k

    sbatch --partition standard --wrap "/usr/bin/time -v -o all_ccs_1k_NM1000_dp0.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/customized/XRSBK_20221007_S64018_PL100268287-1_C01.ccs.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa \
 work_dir=/home/shuaiw/methylation/data/borg/bench/all_ccs_1k_dp0 read_type=ccs min_len=1000 max_NM=1000 min_cov=0"\
  --job-name=all_dp0


    sbatch --partition standard --wrap "/usr/bin/time -v -o all_subreads2.time snakemake \
--config whole_bam=/home/shuaiw/methylation/data/borg/all_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
whole_ref=/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa\
 work_dir=/home/shuaiw/methylation/data/borg/bench/all_subreads2 min_len=1000 max_NM=100000 min_cov=3"\
  --job-name=all

  snakemake --config whole_bam=/home/shuaiw/methylation/data/borg/b_contigs/11.align.bam\
 whole_ref=/home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa\
  work_dir=/home/shuaiw/methylation/data/borg/new_test11 clean=False


  runMetaBat.sh --noBinOut --saveCls all_break.contigs.fa /home/shuaiw/borg/all_break/all_break.align.ccs.bam

