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

# ref_dir=/home/shuaiw/borg/bench/soil/run2/contigs/
# bam_dir=/home/shuaiw/borg/bench/soil/run2/bams/

ref_dir=/home/shuaiw/borg/bench/all_subreads2/contigs/
bam_dir=/home/shuaiw/borg/bench/all_subreads2/bams/
ipd_sum_dir=/home/shuaiw/borg/bench/soil/run2/ipdSum/

params_list=()

for file in $ref_dir/*.fa; do
    # echo "Processing file: $file"
    ## get the basename of the file
    base=$(basename "$file" .fa)
    bam=$bam_dir/$base.bam
    prefix=$ipd_sum_dir/$base
    # echo $bam $file $prefix
    params_list+=("$bam $file $prefix")
done


# # Print the params_list to verify
# for param in "${params_list[@]}"; do
#     echo "$param"
# done


# Run in parallel: each param set is split into three arguments for the script
parallel --colsep ' ' bash /home/shuaiw/Methy/benchmark/compare/run_ipdSummary.sh {1} {2} {3} ::: "${params_list[@]}"
