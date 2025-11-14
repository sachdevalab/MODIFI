while read -r genome_path; do 
  if [ -f "$genome_path" ]; then 
    cp "$genome_path" /home/shuaiw/borg/paper/isolation/GTDB_tree/meta_genomes/
    echo "Copied $(basename "$genome_path")"
  else 
    echo "File not found: $genome_path"
  fi
done < /home/shuaiw/borg/paper/specificity/genome.list
