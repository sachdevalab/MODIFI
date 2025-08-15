#!/bin/bash

# my_prefix=cow_bioreactor_1
# fq1=/groups/diamond/sequences/2023/RuReacBro_20230708_11_72h_200ppm_r1_HiC/raw.d/RuReacBro_20230708_11_72h_200ppm_r1_HiC_trim_clean.PE.1.fastq.gz
# fq2=/groups/diamond/sequences/2023/RuReacBro_20230708_11_72h_200ppm_r1_HiC/raw.d/RuReacBro_20230708_11_72h_200ppm_r1_HiC_trim_clean.PE.2.fastq.gz

threads=40
### cow_bioreactor_1 RuReacBro_20230708_11_72h_200ppm_r1_LR /groups/diamond/sequences/2023/RuReacBro_20230708_11_72h_200ppm_r1_HiC/raw.d/RuReacBro_20230708_11_72h_200ppm_r1_HiC_trim_clean.PE.1.fastq.gz
### cow_bioreactor_2 RuReacBro_20230708_26_72h_NC_r4_LR
### cow_bioreactor_3 RuReacBro_20230708_Comb_RF_HMW_LR
### cow_bioreactor_4  RuReacBro_20230708_9_72h_NC_r2_HMW_LR /groups/diamond/sequences/2023/RuReacBro_20230708_9_72h_NC_r2_HiC/raw.d/RuReacBro_20230708_9_72h_NC_r2_HiC_trim_clean.PE.1.fastq.gz
### cow_bioreactor_5 RuReacBro_20230708_12_72h_200ppm_r2_HMW_LR /groups/diamond/sequences/2023/RuReacBro_20230708_12_72h_200ppm_r2_HiC/raw.d/RuReacBro_20230708_12_72h_200ppm_r2_HiC_trim_clean.PE.1.fastq.gz
### cow_1 RuReacBro_20230708_Comb_RF_LR.

my_prefix=cow_bioreactor_4
fq1=/groups/diamond/sequences/2023/RuReacBro_20230708_9_72h_NC_r2_HiC/raw.d/RuReacBro_20230708_9_72h_NC_r2_HiC_trim_clean.PE.1.fastq.gz
fq2=/groups/diamond/sequences/2023/RuReacBro_20230708_9_72h_NC_r2_HiC/raw.d/RuReacBro_20230708_9_72h_NC_r2_HiC_trim_clean.PE.2.fastq.gz

ref=/home/shuaiw/borg/paper/run2/${my_prefix}/${my_prefix}.hifiasm.p_ctg.rename.fa
workdir=/home/shuaiw/borg/paper/run2/${my_prefix}/hic
prefix=$workdir/${my_prefix}

mkdir -p $workdir
bwa index $ref

bwa mem -t $threads -5SP $ref $fq1 $fq2 | \
samtools view -S -h -b -q 30 -F 2316 -@ $threads| \
samtools sort -n -@ $threads -o ${prefix}_hic.bam

/home/shuaiw/bin/bin3C/bin3C.py mkmap -v  $ref ${prefix}_hic.bam $workdir/bin3c -e MluCI
/home/shuaiw/bin/bin3C/bin3C.py cluster -v $workdir/bin3c/contact_map.p.gz $workdir/bin3c_clust


my_prefix=cow_bioreactor_5
fq1=/groups/diamond/sequences/2023/RuReacBro_20230708_12_72h_200ppm_r2_HiC/raw.d/RuReacBro_20230708_12_72h_200ppm_r2_HiC_trim_clean.PE.1.fastq.gz
fq2=/groups/diamond/sequences/2023/RuReacBro_20230708_12_72h_200ppm_r2_HiC/raw.d/RuReacBro_20230708_12_72h_200ppm_r2_HiC_trim_clean.PE.2.fastq.gz

ref=/home/shuaiw/borg/paper/run2/${my_prefix}/${my_prefix}.hifiasm.p_ctg.rename.fa
workdir=/home/shuaiw/borg/paper/run2/${my_prefix}/hic
prefix=$workdir/${my_prefix}

mkdir -p $workdir
bwa index $ref

bwa mem -t $threads -5SP $ref $fq1 $fq2 | \
samtools view -S -h -b -q 30 -F 2316 -@ $threads| \
samtools sort -n -@ $threads -o ${prefix}_hic.bam

/home/shuaiw/bin/bin3C/bin3C.py mkmap -v  $ref ${prefix}_hic.bam $workdir/bin3c -e MluCI
/home/shuaiw/bin/bin3C/bin3C.py cluster -v $workdir/bin3c/contact_map.p.gz $workdir/bin3c_clust


# bwa mem -t \$SLURM_CPUS_ON_NODE -5SP \$ref \$fq1 \$fq2 | \
# samtools view -S -h -b -q 30 -F 2316 -@ \$SLURM_CPUS_ON_NODE| \
# samtools sort -@ \$SLURM_CPUS_ON_NODE -o ${prefix}_hic.bam

# bwa index $ref

# bwa mem -t $threads -5SP $ref $fq1 $fq2 | \
# samtools view -S -h -b -q 30 -F 2316 -@ $threads| \
# samtools sort -n -@ $threads -o ${prefix}_hic.bam
# # samtools index ${prefix}_hic.bam

# # # Mark duplicate alignments.
# # samtools collate -@ 64 -O -u ${prefix}_hic.bam | samtools fixmate -@ 64 -m -u - - | samtools sort -@ 64 -u - | samtools markdup -@ 64 - ${prefix}_hic_markdup.bam
# # # Sort by name and make BAM text file.
# # samtools sort -@ 64 -n ${prefix}_hic_markdup.bam | samtools view -@ 64 - > ${prefix}_hic_markdup_align.tsv

# /home/shuaiw/bin/bin3C/bin3C.py mkmap -v  $ref ${prefix}_hic.bam $workdir/bin3c -e MluCI
# /home/shuaiw/bin/bin3C/bin3C.py cluster -v $workdir/bin3c/contact_map.p.gz $workdir/bin3c_clust

# coverm genome --bam-files ${prefix}_hic_markdup.bam --genome-fasta-directory /groups/diamond/projects/animal/rumen/RuReacBro20203/HiC/nr_bins_circular_elements/sequence --genome-fasta-extension fa --min-covered-fraction 0 --output-format sparse --methods length count mean variance covered_fraction covered_bases -t 12 --output-file ${prefix}_hic_markdup_coverm_genome.out
# coverm contig --bam-files ${prefix}_hic.bam --min-covered-fraction 0 --output-format sparse --methods length count mean variance covered_fraction covered_bases -t 12 --output-file ${prefix}_hic_markdup_coverm_contig.out
