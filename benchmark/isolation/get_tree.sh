#!/bin/bash

# Create symbolic links for all hifiasm FASTA files
source_dir="/groups/banfield/projects/multienv/methylation_temp/batch2_results"
input_dir="/groups/banfield/projects/multienv/methylation_temp/GTDB_tree/genomes"
out_dir="/groups/banfield/projects/multienv/methylation_temp/GTDB_tree/gtdb_tree2"

# Create the genomes directory if it doesn't exist
mkdir -p $input_dir

# Create symbolic links for all *.hifiasm.p_ctg.rename.fa files
echo "Creating symbolic links..."
for fasta_file in $source_dir/*/*.hifiasm.p_ctg.rename.fa; do
    if [ -f "$fasta_file" ]; then
        # Extract the sample ID from the path (e.g., ERR10357060)
        sample_id=$(basename $(dirname "$fasta_file"))
        
        # Create symbolic link with .fa extension
        target_link="$input_dir/${sample_id}.fa"
        
        if [ ! -e "$target_link" ]; then
            ln -s "$fasta_file" "$target_link"
            echo "Created link: $target_link -> $fasta_file"
        else
            echo "Link already exists: $target_link"
        fi
    fi
done

echo "Symbolic link creation completed."
echo "Found $(ls -1 $input_dir/*.fa 2>/dev/null | wc -l) FASTA files in $input_dir"


gtdbtk  identify --genome_dir $input_dir \
        --out_dir $out_dir \
        --extension fa \
        --cpus 32 \

gtdbtk  align --identify_dir $out_dir \
        --out_dir $out_dir \
        --cpus 32 \

gtdbtk infer --msa_file $out_dir/align/gtdbtk.bac120.user_msa.fasta.gz \
        --out_dir $out_dir --cpus 32

