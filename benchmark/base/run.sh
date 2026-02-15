python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/pure/native \
  --whole_bam /home/shuaiw/methylation/data/published_data/fanggang/C227/native.align.bam \
  --whole_ref /home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa \
  --read_type subreads \
  --min_cov 0 \
  --min_score 0 \
  --kmer_mean_db /home/shuaiw/borg/paper/base/pure/control/control/control_db.up7.down3.mean.dat \
  --kmer_num_db /home/shuaiw/borg/paper/base/pure/control/control/control_db.up7.down3.num.dat \
  --threads 10 

python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/meta/native \
  --whole_bam /home/shuaiw/borg/allison/ecoli/soil_p0.01_C227_align.align.bam \
  --whole_ref /home/shuaiw/borg/allison/ecoli/soil_ecoli.fa \
  --read_type subreads \
  --min_cov 0 \
  --min_score 0 \
  --threads 10

python /home/shuaiw/mGlu/main.py \
  --work_dir /home/shuaiw/borg/paper/base/meta/control \
  --whole_bam /home/shuaiw/borg/allison/ecoli/soil_p0.01_C227_WGA_align.align.bam \
  --whole_ref /home/shuaiw/borg/allison/ecoli/soil_ecoli.fa \
  --read_type subreads \
  --min_cov 0 \
  --min_score 0 \
  --threads 10
