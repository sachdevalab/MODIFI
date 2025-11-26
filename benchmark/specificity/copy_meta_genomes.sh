#!/bin/bash

### shell script copy genome file list in /home/shuaiw/borg/paper/specificity/genome.list
### to the folder /home/shuaiw/borg/paper/isolation/GTDB_tree/meta_genomes/

### genome list file contain on file path per line

# Set variables
GENOME_LIST="/home/shuaiw/borg/paper/specificity/genome.list"
DEST_DIR="/home/shuaiw/borg/paper/isolation/GTDB_tree/meta_genomes/"

# Create destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Check if genome list file exists
if [ ! -f "$GENOME_LIST" ]; then
    echo "Error: Genome list file $GENOME_LIST not found!"
    exit 1
fi

# Counter for tracking progress
count=0
total=$(wc -l < "$GENOME_LIST")

echo "Starting to copy $total genome files..."
echo "Source list: $GENOME_LIST"
echo "Destination: $DEST_DIR"
echo ""

# Read each line (file path) from the genome list
while IFS= read -r genome_path; do
    # Skip empty lines
    if [ -z "$genome_path" ]; then
        continue
    fi
    
    # Check if source file exists
    if [ -f "$genome_path" ]; then
        # Get just the filename from the full path
        filename=$(basename "$genome_path")
        
        # Copy the file to destination
        cp "$genome_path" "$DEST_DIR"
        
        # Check if copy was successful
        if [ $? -eq 0 ]; then
            ((count++))
            echo "[$count/$total] Copied: $filename"
        else
            echo "[$count/$total] ERROR copying: $genome_path"
        fi
    else
        echo "WARNING: File not found: $genome_path"
    fi
    
done < "$GENOME_LIST"

echo ""
echo "Copy operation completed!"
echo "Successfully copied $count out of $total files to $DEST_DIR"


