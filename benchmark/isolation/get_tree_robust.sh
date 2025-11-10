#!/bin/bash

# More robust GTDB-Tk pipeline with problematic genome handling
source_dir="/groups/banfield/projects/multienv/methylation_temp/batch2_results"
input_dir="/home/shuaiw/borg/paper/isolation//GTDB_tree/genomes"
out_dir="/home/shuaiw/borg/paper/isolation//GTDB_tree/gtdb_tree3"

# Create the genomes directory if it doesn't exist
mkdir -p $input_dir

# Known problematic genomes (add more as needed)
problematic_genomes=("ERR13656543")

echo "Creating symbolic links (excluding problematic genomes)..."
total_files=0
excluded_files=0

for fasta_file in $source_dir/*/*.hifiasm.p_ctg.rename.fa; do
    if [ -f "$fasta_file" ] && [ -s "$fasta_file" ]; then
        # Check if file has FASTA content (starts with >)
        if ! head -1 "$fasta_file" | grep -q "^>"; then
            echo "Skipping invalid FASTA file (no header): $fasta_file"
            continue
        fi
        
        # Extract the sample ID from the path
        sample_id=$(basename $(dirname "$fasta_file"))
        
        # Check if this genome is in the problematic list
        skip_genome=false
        for prob_genome in "${problematic_genomes[@]}"; do
            if [ "$sample_id" == "$prob_genome" ]; then
                echo "Skipping problematic genome: $sample_id"
                skip_genome=true
                excluded_files=$((excluded_files + 1))
                break
            fi
        done
        
        if [ "$skip_genome" = false ]; then
            # Create symbolic link with .fa extension
            target_link="$input_dir/${sample_id}.fa"
            
            if [ ! -e "$target_link" ]; then
                ln -s "$fasta_file" "$target_link"
                echo "Created link: $target_link -> $fasta_file"
            else
                echo "Link already exists: $target_link"
            fi
            total_files=$((total_files + 1))
        fi
    fi
done

echo "Symbolic link creation completed."
echo "Total FASTA files: $total_files"
echo "Excluded problematic files: $excluded_files"
echo "Files ready for GTDB-Tk: $(ls -1 $input_dir/*.fa 2>/dev/null | wc -l)"


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