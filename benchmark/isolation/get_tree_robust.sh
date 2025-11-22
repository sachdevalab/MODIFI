#!/bin/bash

input_dir="/home/shuaiw/borg/paper/isolation//GTDB_tree/genomes"
out_dir="/home/shuaiw/borg/paper/isolation//GTDB_tree/gtdb_tree3"




# Step 1: Run GTDB-Tk identify
echo "Running GTDB-Tk identify..."
gtdbtk identify --genome_dir $input_dir \
        --out_dir $out_dir \
        --extension fa \
        --cpus 64

# Check results
identify_dir="$out_dir/identify/intermediate_results/marker_genes"
if [ -d "$identify_dir" ]; then
    successful_genomes=$(find "$identify_dir" -name "*_protein.faa" | wc -l)
    echo "Successfully processed $successful_genomes genomes in identify step"
else
    echo "Error: Identify step failed completely"
    exit 1
fi

# Step 2: Run GTDB-Tk align
echo "Running GTDB-Tk align..."
gtdbtk align --identify_dir $out_dir \
        --out_dir $out_dir \
        --cpus 64

if [ $? -ne 0 ]; then
    echo "Error: GTDB-Tk align failed"
    exit 1
fi

# Step 3: Run GTDB-Tk infer
echo "Running GTDB-Tk infer..."
gtdbtk infer --msa_file $out_dir/align/gtdbtk.bac120.user_msa.fasta.gz \
        --out_dir $out_dir --cpus 64

if [ $? -eq 0 ]; then
    echo "GTDB-Tk pipeline completed successfully!"
    echo "Results available in: $out_dir"
else
    echo "Error: GTDB-Tk infer failed"
    exit 1
fi