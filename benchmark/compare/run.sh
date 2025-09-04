#  sbatch  --partition standard --wrap "/usr/bin/time -v -o ipdSummary.time snakemake -j 64"  --job-name=ipdSum 



  sbatch  --partition standard --wrap "python /home/shuaiw/Methy/main.py \
  --work_dir /home/shuaiw/borg/paper/ipdsummary/ERR12723529_mice_our \
  --whole_bam /home/shuaiw/borg/paper/ipdsummary/ERR12723529_mice//ERR12723529_mice.subreads.align.bam \
  --whole_ref /home/shuaiw/borg/paper/run2/ERR12723529_mice/ERR12723529_mice.hifiasm.p_ctg.rename.fa \
  --read_type subreads \
  --min_len 1000 \
  --min_cov 1 \
  --min_frac 0.3 \
  --min_score 30 \
  --min_sites 30 \
  --min_iden 0.85 \
  --threads 64 --run_steps host \
  --mge_file /home/shuaiw/borg/paper/run2/ERR12723529_mice/all_mge.tsv" \
  --job-name=mice_29


  sbatch  --partition standard --wrap "python /home/shuaiw/Methy/main.py \
  --work_dir /home/shuaiw/borg/paper/ipdsummary/ERR12723528_mice_our \
  --whole_bam /home/shuaiw/borg/paper/ipdsummary/ERR12723528_mice/ERR12723528_mice.subreads.align.bam \
  --whole_ref /home/shuaiw/borg/paper/run2/ERR12723528_mice/ERR12723528_mice.hifiasm.p_ctg.rename.fa \
  --read_type subreads \
  --min_len 1000 \
  --min_cov 1 \
  --min_frac 0.3 \
  --min_score 30 \
  --min_sites 30 \
  --min_iden 0.85 \
  --threads 64 --run_steps host \
  --mge_file /home/shuaiw/borg/paper/run2/ERR12723528_mice/all_mge.tsv" \
  --job-name=mice_28
