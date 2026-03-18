#!/bin/bash
#SBATCH --job-name=ipd_29      # Job name
#SBATCH --partition=standard # Partition name

subreads_bam=/home/shuaiw/borg/paper/aws/subreads/m64079_230801_145304.subreads.demux.bc1008_BAK8A_OA__bc1008_BAK8A_OA.bam
ref=/home/shuaiw/borg/paper/run2/ERR12723529_mice/ERR12723529_mice.hifiasm.p_ctg.rename.fa
threads=$SLURM_CPUS_ON_NODE
outdir=/home/shuaiw/borg/paper/ipdsummary/ERR12723529_mice/

prefix=$outdir/ERR12723529_mice
align_bam=$prefix.subreads.align.bam


## construct output directory if not exists
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

# /usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align $ref $subreads_bam $align_bam \
#  --sort -j $threads -J $threads
# samtools index $align_bam
# /home/shuaiw//smrtlink/pbindex $align_bam

/usr/bin/time -v -o $prefix.pbmm2.align.time ~/smrtlink/pbmm2 align --preset SUBREAD -j $SLURM_CPUS_ON_NODE $ref $subreads_bam $prefix.raw.bam 
samtools sort -T $prefix -@ $SLURM_CPUS_ON_NODE -o $align_bam $prefix.raw.bam
rm $prefix.raw.bam
samtools index $align_bam
/home/shuaiw//smrtlink/pbindex $align_bam


/usr/bin/time -v -o $prefix.ipdSummary.time \
~/smrtlink/ipdSummary $align_bam --reference $ref --debug --numWorkers $threads \
 --gff $prefix.gff --csv $prefix.csv  --methylFraction --outfile $prefix


/usr/bin/time -v -o $prefix.motifmaker.time ~/smrtlink/motifMaker find -f $ref -g $prefix.gff -j $threads -o $prefix.motif.csv -m 30

/usr/bin/time -v -o $prefix.reprocess.time ~/smrtlink/motifMaker reprocess -m $prefix.motif.csv \
 -f $ref -g $prefix.gff -c $prefix.csv -o $prefix.reprocess.gff