### subreads to hifi


hifi_bam=/home/shuaiw/borg/m84039_230624_013044_s3.hifi_reads.bc2021.bam
prefix=soil_s1_2021
work_dir=/home/shuaiw/borg/assembly/soil/$prefix

## build the dir
mkdir -p $work_dir

## bam to fastq using samtools
samtools fastq $hifi_bam > $work_dir/${prefix}.hifi.fastq
## gzip the fastq
gzip $work_dir/${prefix}.hifi.fastq

## assembly using hifiasm
# hifiasm_meta -o $work_dir/${prefix}.hifiasm -t 64  $work_dir/${prefix}.hifi.fastq.gz  > $work_dir/${prefix}.hifiasm.log
## transform to fasta


## identify MGE 

## map reads to assembly

