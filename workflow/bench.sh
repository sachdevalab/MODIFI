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
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_NM3 \
  read_type=ccs min_len=1000 max_NM=3"\
   --job-name=zymo2

     sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.10pct.bam \
 whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.1 \
  read_type=ccs min_len=1000 max_NM=3 min_cov=1"\
   --job-name=zymo2


    sbatch --partition standard --wrap " genomad end-to-end --relaxed --cleanup --enable-score-calibration \
  --threads 64 --sensitivity 7.0 --force-auto \
  /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_genomad/ \
  /groups/diamond/databases/genomad/v1.7/" --job-name=genomad


  sbatch  --partition standard --wrap "snakemake --config \
 whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.align.bam \
 whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo \
  read_type=ccs clean=False"\
   --job-name=zymo2

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1 \
  read_type=ccs min_len=1000 max_NM=3 min_cov=1 clean=False min_frac=0.1 min_score=13" \
  --job-name=zymo1

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov5 \
  read_type=ccs min_len=1000 max_NM=3 min_cov=5 clean=False min_frac=0.1 min_score=13" \
  --job-name=zymo2

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov5_s30 \
  read_type=ccs min_len=1000 max_NM=3 min_cov=5 clean=False min_frac=0.1 min_score=30" \
  --job-name=zymo3

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s20 \
  read_type=ccs min_len=1000 max_NM=3 min_cov=1 clean=False min_frac=0.1 min_score=20" \
  --job-name=zymo4

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s10 \
  read_type=ccs min_len=1000 max_NM=3 min_cov=1 clean=False min_frac=0.1 min_score=30 min_sites=10" \
  --job-name=zymo4

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30 \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 clean=False min_frac=0.1 min_score=30 min_sites=10" \
  --job-name=zymo4

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_gmm \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 clean=False min_frac=0.1 min_score=30 min_sites=10" \
  --job-name=zymo4


sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.10pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.1_cov1_s30_filter2 \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 clean=False" \
  --job-name=zymo4

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_filter \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10" \
  --job-name=zymo4

sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_filter2 \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 clean=False" \
  --job-name=p5


### try segment with dp 5
 sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_rec3 \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 clean=False \
  plasmid_file=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list" \
  --job-name=recover    

## test set new motif critera's impact on the result
 sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_rec3 \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=50 clean=False \
  plasmid_file=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list" \
  --job-name=recover    


sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/merged2_break.align.ccs.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_break.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_break2 \
  read_type=ccs min_len=1000 max_NM=3" \
  --job-name=zymo_break



sbatch  --partition standard --wrap "hifiasm_meta -o /home/shuaiw/methylation/data/borg/bench/zymop5_ass/zymop5 \
-t 64 /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.fq.gz"\
  --job-name=assemble

samtools fastq /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam >/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.fq 


sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20.time snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_200ppm_r1_LR_scaffold.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_new \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10"\
  --job-name=pf



python comp_ipd_ratio.py /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/ipd_ratio/E_coli_H10407_1.ipd3.csv \
  /home/shuaiw/borg/bench/test/test.csv \
  /home/shuaiw/borg/bench/test/E_coli_H10407_1.gff \
  /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/contigs/E_coli_H10407_1.fa \
  /home/shuaiw/borg/bench/test/test.png \
  1

/home/shuaiw/smrtlink/motifMaker find -g /home/shuaiw/borg/bench/test/E_coli_H10407_1.gff \
      -o /home/shuaiw/borg/bench/test/test.motif \
          -m 30 -j 10 -f /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/contigs/E_coli_H10407_1.fa

checkm2 predict --input  bins/ --output-directory  test --force -x .fasta


sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20.time snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_200ppm_r1_LR_scaffold.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_new_long \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 clean=False \
  plasmid_file=/groups/diamond/projects/animal/rumen/RuReacBro20203/assembly/RuReacBro_20230708_11_72h_200ppm_r1_LR/genomad/circular_elements_summary/circular_elements_plasmid_summary.tsv \
  --rerun-incomplete"\
  --job-name=pf

 sbatch  --partition standard --wrap " snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/align/RuReacBro_20230708_12_72h_200ppm_r2_HMW_LR.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_12_72h_200ppm_r2_HMW_LR_bin \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 \
  plasmid_file=/home/shuaiw/borg/pengfan/contigs/MGE.list"\
  --job-name=pf3

 sbatch  --partition standard --wrap " snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/align/RuReacBro_20230708_26_72h_NC_r4_LR.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_26_72h_NC_r4_LR_bin \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 \
  plasmid_file=/home/shuaiw/borg/pengfan/contigs/MGE.list"\
  --job-name=pf4

 sbatch  --partition standard --wrap " snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/align/RuReacBro_20230708_9_72h_NC_r2_HMW_LR.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_9_72h_NC_r2_HMW_LR_bin \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 \
  plasmid_file=/home/shuaiw/borg/pengfan/contigs/MGE.list"\
  --job-name=pf5

 sbatch  --partition standard --wrap " snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/align/RuReacBro_20230708_Comb_RF_HMW_LR.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_Comb_RF_HMW_LR_bin \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 \
  plasmid_file=/home/shuaiw/borg/pengfan/contigs/MGE.list"\
  --job-name=pf6

 sbatch  --partition standard --wrap " snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/align/RuReacBro_20230708_Comb_RF_LR.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_Comb_RF_LR_bin \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 \
  plasmid_file=/home/shuaiw/borg/pengfan/contigs/MGE.list"\
  --job-name=pf7


python /home/shuaiw/Methy/motif_profile.py /home/shuaiw/methylation/data/borg/bench/all_ccs_new/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_435_C.fa /home/shuaiw/methylation/data/borg/bench/all_ccs_new/gffs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_435_C.gff /home/shuaiw/methylation/data/borg/bench/all_ccs_new/all.motifs.csv /home/shuaiw/methylation/data/borg/bench/all_ccs_new/profiles/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_435_C.motifs.profile.csv /home/shuaiw/methylation/data/borg/bench/all_ccs_new/ipd_ratio/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_435_C.ipd3.csv  0.1 10 30 1

## test set new motif critera's impact on the result
 sbatch  --partition standard --wrap "snakemake --config \
  whole_bam=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  whole_ref=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  work_dir=/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_rec3 \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=50 clean=False \
  plasmid_file=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list" \
  --job-name=recover    

sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/methylation/data/borg/bench/pipeline_zymo.time python main.py \
  --work_dir /home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_rec5 \
  --whole_bam /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list  " \
  --job-name=pipeline

   sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20.new.time snakemake --config \
  whole_bam=/home/shuaiw/borg/pengfan/align/RuReacBro_20230708_11_72h_20_new.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  work_dir=/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.1 min_score=30 min_sites=10 clean=False \
  plasmid_file=/home/shuaiw/borg/pengfan/contigs/MGE.list"\
  --job-name=pf2

python main.py \
  --work_dir /home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin \
  --whole_bam /home/shuaiw/borg/pengfan/align/RuReacBro_20230708_11_72h_20_new.align.ccs.bam \
  --whole_ref /home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid_file /home/shuaiw/borg/pengfan/contigs/MGE.list \
  --bin_file  /home/shuaiw/borg/pengfan/10mgs_bins.tab \
  --run_steps merge \
  --threads 20 

python cal_invasion_score.py \
  --work_dir /home/shuaiw/methylation/data/borg/pengfan/RuReacBro_20230708_11_72h_20_bin\
  --bin_file  /home/shuaiw/borg/pengfan/10mgs_bins.tab\
  --min_frac 0.4 \
  --threads 10 \
  --plasmid_file /home/shuaiw/borg/pengfan/contigs/MGE.list

python cal_invasion_score.py \
  --work_dir /home/shuaiw/methylation/data/borg/pengfan/RuReacBro_20230708_11_72h_20_bin\
  --bin_file  /home/shuaiw/borg/pengfan/10mgs_bins.tab\
  --min_frac 0.4 --threads 10 --plasmid RuReacBro_20230708_26_72h_NC_r4_LR_9420_C

 sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/bench/soil_zymo/time snakemake --config \
  whole_bam=/home/shuaiw/borg/bench/soil_zymo/soil_zymo.align.ccs.bam \
  whole_ref=/home/shuaiw/borg/contigs/soil_zymo.fa \
  work_dir=/home/shuaiw/borg/bench/soil_zymo/run \
  read_type=ccs min_len=1000 max_NM=3000 min_cov=1 min_frac=0.4 min_score=30 min_sites=30 clean=False \
  plasmid_file=/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list"\
  --job-name=soil

 sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/bench/soil_zymo/run/run.time python main.py \
  --work_dir /home/shuaiw/borg/bench/soil_zymo/run \
  --whole_bam /home/shuaiw/borg/bench/soil_zymo/soil_zymo.align.ccs.bam \
  --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list \
  --run_steps profile merge host \
  --threads 64 "\
  --job-name=soil

 sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/borg/bench/soil_zymo/run3.time python main.py \
  --work_dir /home/shuaiw/borg/bench/soil_zymo/run3 \
  --whole_bam /home/shuaiw/borg/bench/soil_zymo/soil_zymo.align.ccs.bam \
  --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 10 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list \
  --threads 64 "\
  --job-name=soil3


 sbatch  --partition standard --wrap "python main.py \
  --work_dir /home/shuaiw/borg/bench/soil_zymo/run2 \
  --whole_bam /home/shuaiw/borg/bench/soil_zymo/soil_zymo.align.ccs.bam \
  --whole_ref /home/shuaiw/borg/contigs/soil_zymo.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list \
  --run_steps profile merge host\
  --threads 64 "\
  --job-name=soil2


sbatch  --partition standard --wrap "python main.py \
  --work_dir /home/shuaiw/borg/allison/ecoli/native \
  --whole_bam /home/shuaiw/borg/allison/ecoli/soil_p0.01_C227_align.align.bam \
  --whole_ref /home/shuaiw/borg/allison/ecoli/soil_ecoli.fa \
  --read_type subreads \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --run_steps control compare motif profile merge \
  --threads 64 "\
  --job-name=c227

sbatch  --partition standard --wrap "python main.py \
  --work_dir /home/shuaiw/borg/allison/ecoli/WGA \
  --whole_bam /home/shuaiw/borg/allison/ecoli/soil_p0.01_C227_WGA_align.align.bam \
  --whole_ref /home/shuaiw/borg/allison/ecoli/soil_ecoli.fa \
  --read_type subreads \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --run_steps profile merge \
  --threads 64 "\
  --job-name=c227_WGA


sbatch  --partition standard --wrap " python main.py \
  --work_dir /home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2 \
  --whole_bam /home/shuaiw/borg/pengfan/align/RuReacBro_20230708_11_72h_20_new.align.ccs.bam \
  --whole_ref /home/shuaiw/borg/pengfan/contigs/nr_bins_circular_elements.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid_file /home/shuaiw/borg/pengfan/contigs/MGE.list \
  --bin_file  /home/shuaiw/borg/pengfan/10mgs_bins.tab \
  --run_steps profile merge \
  --threads 64"\
  --job-name=new_pf

sbatch  --partition standard --wrap " python main.py \
  --work_dir /home/shuaiw/methylation/data/borg/bench/mock3 \
  --whole_bam /home/shuaiw/methylation/data/published_data/fanggang/align/Mock_JF8.align.bam \
  --whole_ref /home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa \
  --read_type subreads \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid PDZQ01000142.1\
  --threads 64"\
  --job-name=J8_mock



sbatch  --partition standard --wrap "/usr/bin/time -v -o /home/shuaiw/methylation/data/borg/bench/pipeline_zymo.time python main.py \
  --work_dir /home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30_rec7 \
  --whole_bam /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/align/m64004_210929_143746.5pct.bam \
  --whole_ref /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list\
  --threads 10  " \
  --job-name=pipeline

