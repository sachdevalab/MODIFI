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

# Step 1: Run GTDB-Tk identify
echo "Running GTDB-Tk identify..."
gtdbtk identify --genome_dir $input_dir \
        --out_dir $out_dir \
        --extension fa \
        --cpus 32

# Check if identify step completed successfully
if [ $? -ne 0 ]; then
    echo "Error: GTDB-Tk identify failed. Check the logs above."
    echo "Some genomes may have failed during Prodigal analysis."
    echo "Continuing with available genomes..."
fi

# Check what genomes were successfully processed
identify_dir="$out_dir/identify/intermediate_results/marker_genes"
if [ -d "$identify_dir" ]; then
    successful_genomes=$(find "$identify_dir" -name "*_protein.faa" | wc -l)
    echo "Successfully processed $successful_genomes genomes in identify step"
    
    # List failed genomes for debugging
    echo "Checking for failed genomes..."
    for genome_dir in "$identify_dir"/*/; do
        if [ -d "$genome_dir" ]; then
            genome_name=$(basename "$genome_dir")
            protein_file="${genome_dir}/${genome_name}_protein.faa"
            if [ ! -f "$protein_file" ]; then
                echo "Failed genome: $genome_name (missing protein file)"
            fi
        fi
    done
else
    echo "Error: Identify directory not found at $identify_dir"
    exit 1
fi

# Step 2: Run GTDB-Tk align (only if identify completed)
echo "Running GTDB-Tk align..."
gtdbtk align --identify_dir $out_dir \
        --out_dir $out_dir \
        --cpus 32

# Check if align step completed successfully
if [ $? -ne 0 ]; then
    echo "Error: GTDB-Tk align failed. Check the logs above."
    exit 1
fi

# Step 3: Run GTDB-Tk infer (only if align completed)
echo "Running GTDB-Tk infer..."
gtdbtk infer --msa_file $out_dir/align/gtdbtk.bac120.user_msa.fasta.gz \
        --out_dir $out_dir --cpus 32

if [ $? -eq 0 ]; then
    echo "GTDB-Tk pipeline completed successfully!"
else
    echo "Error: GTDB-Tk infer failed. Check the logs above."
    exit 1
fi

