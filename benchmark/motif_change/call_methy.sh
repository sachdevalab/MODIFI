outdir=/home/shuaiw/borg/paper/motif_change/new_alignment/
ref=/home/shuaiw/borg/paper/run2/infant_2/infant_2_methylation4/contigs/infant_2_3_C.fa

if [ ! -d $outdir ]; then
    mkdir $outdir
fi


ccs_bam=/home/shuaiw/borg/paper/aws/infant/NANO_2_INF1240040_2G1_pacbio.bam

for week in week6 week5 week4 week3; do
    prefix=$outdir/$week
    align_bam=$prefix.align.ccs.bam

    python /home/shuaiw/mGlu/main.py \
      --work_dir $prefix \
      --whole_bam $align_bam \
      --whole_ref $ref \
      --kmer_mean_db /home/shuaiw/methylation/data/borg/paper/run2/infant_2/infant_2_methylation4/control/control_db.up7.down3.mean.dat \
      --kmer_num_db /home/shuaiw/methylation/data/borg/paper/run2/infant_2/infant_2_methylation4/control/control_db.up7.down3.num.dat \
      --threads 10
done