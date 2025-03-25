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
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo4_NM0 \
  read_type=ccs min_len=5000 max_NM=0"\
   --job-name=zymo2


   sbatch  --partition standard --wrap " /home/shuaiw/bin/bin/jasmine -j 64 /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.align.bam \
   /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.align.jasmine.bam"\
   --job-name=jasmine

   sbatch  --partition standard --wrap "/home/shuaiw/bin/pb-CpG-tools-v3.0.0-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores \
    --bam /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.align.jasmine.bam \
     --output-prefix /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746_cpgtools --threads 64"\
   --job-name=cpgtools


   sbatch --partition standard --wrap "snakemake \
--config whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/merge.align.bam \
whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/three_species.fa \
work_dir=/home/shuaiw/methylation/data/borg/bench/merge2" --job-name=merge

   sbatch --partition standard --wrap "snakemake \
--config whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/align/merge_WGA.align.bam \
whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/three_species.fa \
work_dir=/home/shuaiw/methylation/data/borg/bench/merge_WGA" --job-name=merge


      sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/borg/break_contigs/XRSBK_20221007_S64018_PL100268287-1_C01.align.bam \
 whole_ref=/home/shuaiw/borg/contigs/break_contigs.fasta \
  work_dir=/home/shuaiw/methylation/data/borg/bench/breakdb" --job-name=breakdb


        sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/C227/native.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/C227/native\
  kmer_mean_db=/home/shuaiw/borg/bench/breakdb/control/control_db.up7.down3.mean.dat\
  kmer_num_db=/home/shuaiw/borg/bench/breakdb/control/control_db.up7.down3.num.dat" --job-name=native


          sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/C227/WGA.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/C227/WGA2" --job-name=WGA


          sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/published_data/fanggang/C227/native.align.bam \
 whole_ref=/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/C227/native2\
  kmer_mean_db=/home/shuaiw/borg/bench/C227/WGA2/control/control_db.up7.down3.mean.dat\
  kmer_num_db=/home/shuaiw/borg/bench/C227/WGA2/control/control_db.up7.down3.num.dat" --job-name=native


   sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.align.bam \
 whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref \
  read_type=ccs min_len=2000 max_NM=5 min_cov=5"\
   --job-name=zymo2

     sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.align.bam \
 whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo6_NM3 \
  read_type=ccs min_len=1000 max_NM=3"
   --job-name=zymo2

     sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20.time snakemake --config \
 whole_bam=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20.align.ccs.bam \
 whole_ref=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_200ppm_r1_LR_scaffold.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_more \
  read_type=ccs min_len=5000 max_NM=1000 min_cov=3"\
   --job-name=pf

       sbatch  --partition standard --wrap "/usr/bin/time -v -o zymo_linear.time snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/relinear_ref/relinear_ref.align.ccs.bam \
 whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/relinear_ref/merged.fasta \
  work_dir=/home/shuaiw/borg/bench/zymo_linear \
  read_type=ccs min_len=2000 max_NM=100 min_cov=5"\
   --job-name=zymo
  
   sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.align.ccs.bam \
 whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref \
  read_type=ccs min_len=2000 max_NM=5 min_cov=5"\
   --job-name=zymo2

    sbatch --partition standard --wrap " genomad end-to-end --relaxed --cleanup --enable-score-calibration \
  --threads 64 --sensitivity 7.0 --force-auto \
  /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_genomad/ \
  /groups/diamond/databases/genomad/v1.7/" --job-name=genomad
