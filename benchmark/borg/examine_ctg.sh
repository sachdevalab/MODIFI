fq1=/groups/banfield/sequences/2022/SR-VP_9_9_2021_81_5A_0_75m_2/raw.d/SR-VP_9_9_2021_81_5A_0_75m_2_trim_clean.PE.1.fastq.gz
fq2=/groups/banfield/sequences/2022/SR-VP_9_9_2021_81_5A_0_75m_2/raw.d/SR-VP_9_9_2021_81_5A_0_75m_2_trim_clean.PE.2.fastq.gz
ref=/home/shuaiw/borg/paper/circos/borg2/soil_1_11304_C.fa
outdir=/home/shuaiw/borg/paper/borg_data/circular/
prefix=soil_1_11304_C


# bbmap.sh pigz=t unpigz=t ambiguous=random minid=1.96 idfilter=0.97\
#  threads=64 out=stdout.sam editfilter=5 \
#  out=$outdir/$prefix.sam in1=$fq1\
#   in2=$fq2\
#     ref=$ref\
#       nodisk | shrinksam | sambam > $outdir/$prefix.bam

## sort the $outdir/$prefix.sam and output to $outdir/$prefix.bam and index the bam file
## only keep reads with map q >=20
samtools view -b -q 20 -o $outdir/$prefix.bam $outdir/$prefix.sam
samtools sort -o $outdir/$prefix.sorted.bam $outdir/$prefix.bam
samtools index $outdir/$prefix.sorted.bam