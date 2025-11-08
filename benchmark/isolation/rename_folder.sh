#!/bin/bash

# Script to rename *_methylation/ directories to *_methylation2/
# Path: /home/shuaiw/borg/paper/isolation/bacteria/*/*_methylation/

base_dir="/home/shuaiw/borg/paper/isolation/bacteria"

echo "Starting to rename methylation directories..."
echo "Base directory: $base_dir"

# Counter for tracking
count=0
renamed=0
errors=0

# Find all directories matching the pattern *_methylation
for methylation_dir in "$base_dir"/*/*_methylation/; do
    # Check if the directory actually exists (in case no matches found)
    if [ -d "$methylation_dir" ]; then
        count=$((count + 1))
        
        # Remove trailing slash for cleaner processing
        methylation_dir="${methylation_dir%/}"
        
        # Create the new directory name by replacing _methylation with _methylation2
        new_dir="${methylation_dir%_methylation}_methylation2"
        
        echo "Processing: $methylation_dir"
        
        # Check if target directory already exists
        if [ -d "$new_dir" ]; then
            echo "  [WARNING] Target already exists: $new_dir"
            errors=$((errors + 1))
        else
            # Perform the rename
            if mv "$methylation_dir" "$new_dir"; then
                echo "  [SUCCESS] Renamed to: $new_dir"
                renamed=$((renamed + 1))
            else
                echo "  [ERROR] Failed to rename: $methylation_dir"
                errors=$((errors + 1))
            fi
            echo mv "$methylation_dir" "$new_dir";
        fi
    fi
done

# Summary
echo ""
echo "=== SUMMARY ==="
echo "Total directories found: $count"
echo "Successfully renamed: $renamed"
echo "Errors/warnings: $errors"

if [ $count -eq 0 ]; then
    echo "No directories matching pattern '*_methylation/' found in $base_dir"
fi