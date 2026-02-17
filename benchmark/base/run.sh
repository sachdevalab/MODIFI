python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/pure/control \
  --whole_bam /home/shuaiw/methylation/data/published_data/fanggang/C227/WGA.align.bam \
  --whole_ref /home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  --read_type subreads \
  --min_cov 1 \
  --min_score 30 \
  --threads 10

python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/pure/native \
  --whole_bam /home/shuaiw/methylation/data/published_data/fanggang/C227/native.align.bam \
  --whole_ref /home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  --read_type subreads \
  --min_cov 1 \
  --min_score 30 \
  --kmer_mean_db /home/shuaiw/borg/paper/base/pure/control/control/control_db.up7.down3.mean.dat \
  --kmer_num_db /home/shuaiw/borg/paper/base/pure/control/control/control_db.up7.down3.num.dat \
  --threads 10 

python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/meta/native \
  --whole_bam /home/shuaiw/borg/allison/ecoli/soil_p0.01_C227_align.align.bam \
  --whole_ref /home/shuaiw/borg/allison/ecoli/soil_ecoli.fa \
  --read_type subreads \
  --min_cov 1 \
  --min_score 30 \
  --threads 30

python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/meta/control \
  --whole_bam /home/shuaiw/borg/allison/ecoli/soil_p0.01_C227_WGA_align.align.bam \
  --whole_ref /home/shuaiw/borg/allison/ecoli/soil_ecoli.fa \
  --read_type subreads \
  --min_cov 1 \
  --min_score 30 \
  --threads 30

python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/mock3 \
  --whole_bam /home/shuaiw/methylation/data/published_data/fanggang/align/Mock_JF8.align.bam \
  --whole_ref /home/shuaiw/methylation/data/published_data/fanggang/bam/Mock_JF8.fa \
  --read_type subreads \
  --min_len 1000 \
  --min_cov 5 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 100 \
  --threads 32
