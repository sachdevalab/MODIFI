sbatch --partition standard --wrap "snakemake \ 
--config whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/Escherichia_coli_O104:H4_str._C227-11__250_bp_library__native_DNA_.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/bam/GCA_022869985.1_ASM2286998v1_genomic.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/ecoli_native" --job-name=ecoli_native
sbatch --partition standard --wrap "snakemake --config whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/Escherichia_coli_O104:H4_str._C227-11__250_bp_library__whole-genome_amplified_DNA_.align.bam whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/bam/GCA_022869985.1_ASM2286998v1_genomic.fa work_dir=/home/shuaiw/methylation/data/borg/bench/ecoli_control" --job-name=ecoli_control
sbatch --partition standard --wrap "snakemake --config whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/Chromohalobacter_salexigens_DSM-3043__250_bp_library__whole-genome_amplified_DNA_.align.bam whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/bam/GCA_000055785.1_ASM5578v1_genomic.fa work_dir=/home/shuaiw/methylation/data/borg/bench/chr_control" --job-name=chr_control
sbatch --partition standard --wrap "snakemake  \ 
--config whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/Chromohalobacter_salexigens_DSM-3043__250_bp_library__native_DNA_.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/bam/GCA_000055785.1_ASM5578v1_genomic.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/chr_native" --job-name=chr_native
sbatch --partition standard --wrap "snakemake --config whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/merge.align.bam whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/bam/merge.fa work_dir=/home/shuaiw/methylation/data/borg/bench/merge" --job-name=merge


sbatch  --partition standard --wrap "snakemake -s pipeline2.smk --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/Escherichia_coli_O104:H4_str._C227-11__250_bp_library__native_DNA_.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/bam/GCA_022869985.1_ASM2286998v1_genomic.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/ecoli_native2" --job-name=ecoli_native


sbatch  --partition standard --wrap "snakemake -s pipeline2.smk --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/ref/C227_native.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/C227_native" --job-name=ecoli_native

  sbatch  --partition standard --wrap "snakemake -s pipeline2.smk --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/ref/C227_WGA.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/C227_WGA" --job-name=ecoli_WGA

  sbatch  --partition standard --wrap "snakemake -s pipeline2.smk --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/ref/J99_native.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/J99.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/J99_native" --job-name=J99_native

    sbatch  --partition standard --wrap "snakemake -s pipeline2.smk --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/ref/J99_WGA.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/J99.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/J99_WGA" --job-name=J99_WGA


      sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/borg/break_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
 whole_ref=/home/shuaiw/borg/contigs/break_contigs.fasta \
  work_dir=/home/shuaiw/methylation/data/borg/bench/break3" --job-name=break3

        sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/Mock_JF8.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/mock" --job-name=mock

        sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/merge.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/three_species.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/three" --job-name=three

          sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/borg/all_borg2/all_borg2.align.bam \
 whole_ref=/home/shuaiw/methylation/data/borg/all_borg2/all_borgs2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/borg" --job-name=borg

        sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/borg/break_contigs2/XRSBK_20221007_S64018_PL100268287-1_C01.align.ccs.bam \
 whole_ref=/home/shuaiw/borg/contigs/break_contigs.fasta \
  work_dir=/home/shuaiw/methylation/data/borg/bench/break5 max_NM=0 \
  read_type=ccs"\
   --job-name=break4


  sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.align.bam \
 whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo2 \
  read_type=ccs min_len=5000"\
   --job-name=zymo2