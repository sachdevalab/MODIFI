#!/bin/bash

# Define a function that processes two parameters
process_parameters() {
    param1=$1
    param2=$2
    param3=$3
    echo "Processing parameters: $param1 and $param2"
    # Add your processing logic here
}

# Export the function to make it available to parallel
export -f process_parameters

# Define the range of values for generating parameters
start=1
end=5

# Generate the params_list automatically
# params_list=()
# for ((i=start; i<=end; i++)); do
#     params_list+=("param1_value$i param2_value$i")
# done


my_dir=/home/shuaiw/methylation/data/borg/split_bam_dir3/

params_list=()
## for each file in the directory with suffix *fasta
## folder is split_bam_dir2
for file in $my_dir/*.fasta; do
    ## get the basename of the file
    base=$(basename "$file" .fasta)
    bam=$my_dir/$base.bam
    contig=$base
    # echo $bam $file $contig
    params_list+=("$bam $file $contig")
done


# Print the params_list to verify
# for param in "${params_list[@]}"; do
#     echo "$param"
# done


# # Use parallel to process each set of parameters
parallel python split_bam_single.py ::: "${params_list[@]}"
