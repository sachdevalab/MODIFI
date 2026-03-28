
## aligned bam : --whole_bam /home/shuaiw/borg/paper/run2/96plex/96plex.align.bam
## for test
sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
  -o /home/shuaiw/borg/paper/software_test/test_0327/96plex_test_new \
  -b /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/m64004_210929_143746.hifi_reads.bam \
  -r /home/shuaiw/borg/paper/run2/96plex/96plex.hifiasm.p_ctg.rename.fa \
  --read_type hifi \
  --min_len 1000 \
  --min_cov 5 \
  --min_iden 0.97 \
  --min_frac 0.3 \
  --min_score 30 \
  --min_sites 100 \
  --threads 64 \
  --mge_file /home/shuaiw/borg/paper/run2/96plex/all_mge.tsv " \
  --job-name=96plex_test


  sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
  -o /home/shuaiw/borg/paper/software_test/test_0327/96plex_aligned \
  --aligned_bam /home/shuaiw/borg/paper/run2/96plex/96plex.align.bam \
  -r /home/shuaiw/borg/paper/run2/96plex/96plex.hifiasm.p_ctg.rename.fa \
  --read_type hifi \
  --min_len 1000 \
  --min_cov 5 \
  --min_iden 0.97 \
  --min_frac 0.3 \
  --min_score 30 \
  --min_sites 100 \
  --threads 64 \
  --mge_file /home/shuaiw/borg/paper/run2/96plex/all_mge.tsv " \
  --job-name=96plex_test


  sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
  -o /home/shuaiw/borg/paper/software_test/test_0327/ERR12723528_mice_new \
  -b /home/shuaiw/borg/paper/aws/ERR12723528/ERR12723528.ccs.bam \
  -r /home/shuaiw/borg/paper/run2/ERR12723528_mice/ERR12723528_mice.hifiasm.p_ctg.rename.fa \
  --read_type hifi \
  --min_len 1000 \
  --min_cov 5 \
  --min_iden 0.97 \
  --min_frac 0.3 \
  --min_score 30 \
  --min_sites 100 \
  --threads 64 \
  --mge_file /home/shuaiw/borg/paper/run2/ERR12723528_mice/all_mge.tsv " \
  --job-name=ERR12723528_mice

  sbatch  --partition standard --wrap "python /home/shuaiw/MODIFI/main.py \
  -o /home/shuaiw/borg/paper/software_test/test_0327/ERR12723528_mice_aligned \
  --aligned_bam /home/shuaiw/borg/paper/run2/ERR12723528_mice/ERR12723528_mice.align.bam \
  -r /home/shuaiw/borg/paper/run2/ERR12723528_mice/ERR12723528_mice.hifiasm.p_ctg.rename.fa \
  --read_type hifi \
  --min_len 1000 \
  --min_cov 5 \
  --min_iden 0.97 \
  --min_frac 0.3 \
  --min_score 30 \
  --min_sites 100 \
  --threads 64 \
  --mge_file /home/shuaiw/borg/paper/run2/ERR12723528_mice/all_mge.tsv " \
  --job-name=ERR12723528_mice
