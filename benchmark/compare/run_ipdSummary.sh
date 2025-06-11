align_bam=$1
ref=$2
prefix=$3
echo "~/smrtlink/ipdSummary $align_bam --reference $ref --debug --numWorkers 1  --gff $prefix.gff --csv $prefix.csv  --methylFraction --outfile $prefix.out"
~/smrtlink/ipdSummary $align_bam --reference $ref --debug --numWorkers 1  --gff $prefix.gff --csv $prefix.csv  --methylFraction --outfile $prefix.out
