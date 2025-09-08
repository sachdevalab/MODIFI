align_bam=$1
ref=$2
prefix=$3
# echo "~/smrtlink/ipdSummary $align_bam --reference $ref --debug --numWorkers 1  --gff $prefix.gff --csv $prefix.csv  --methylFraction --outfile $prefix.out"
# ~/smrtlink/ipdSummary $align_bam --reference $ref --debug --numWorkers 1  --gff $prefix.gff --csv $prefix.csv  --methylFraction --outfile $prefix.out

~/smrtlink/motifMaker find -f $ref -g $prefix.gff -m 30 -j 10 -o $prefix.motif.csv
~/smrtlink/motifMaker reprocess -c $prefix.csv -f $ref -g $prefix.gff -m $prefix.motif.csv -o $prefix.motif.reprocess.gff