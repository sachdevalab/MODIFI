shopt -s nullglob
for f in /home/shuaiw/borg/paper/run2/*/*_methylation3/host_summary.csv; do
  cp -p "$f" "$(dirname "$f")/host_summary_0911.csv"
done

