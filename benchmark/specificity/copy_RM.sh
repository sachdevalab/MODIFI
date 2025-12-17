#!/bin/bash

# Copy RM_systems directories from methylation3 to methylation4
# For all samples in /home/shuaiw/borg/paper/run2/

base_dir="/home/shuaiw/borg/paper/run2"

for source_dir in ${base_dir}/*/*_methylation3/RM_systems; do
    if [ -d "$source_dir" ]; then
        # Extract the sample directory path
        sample_dir=$(dirname $(dirname "$source_dir"))
        sample_name=$(basename "$sample_dir")
        
        # Construct destination path
        dest_dir="${sample_dir}/${sample_name}_methylation4/RM_systems"
        
        # Copy the RM_systems directory
        echo "Copying: $source_dir -> $dest_dir"
        cp -r "$source_dir" "$dest_dir"
        # break
    fi
done

echo "Done copying RM_systems directories"
