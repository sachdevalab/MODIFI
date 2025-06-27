### subreads to hifi


# hifi_bam=/home/shuaiw/borg/m84039_230624_013044_s3.hifi_reads.bc2021.bam
# prefix=soil_s1_2021
# work_dir=/home/shuaiw/borg/assembly/soil/$prefix

hifi_bam=/home/shuaiw/borg/paper/linkage/m64004_210929_143746.p05.bam
prefix=96plex_p5
work_dir=/home/shuaiw/borg/assembly/96plex/$prefix
threads=20

## build the dir
mkdir -p $work_dir

## bam to fastq using samtools
# samtools fastq $hifi_bam > $work_dir/${prefix}.hifi.fastq
# ## gzip the fastq
# gzip $work_dir/${prefix}.hifi.fastq

## trim reads
# parallel  --plus --eta -j4  “/shared/software/bbmap/39.01/bbduk.sh in={} out={..}.qc.fq.gz minavgquality=20 qtrim=rl trimq=20 threads=5” ::: /raw.d/.fq.gz

## remove host reads, minimap2 map to the reference, collect the reference of hosts, send to rohan , ask allison, add human

# ## assembly using hifiasm
# hifiasm_meta -o $work_dir/${prefix}.hifiasm -t $threads  $work_dir/${prefix}.hifi.fastq.gz  > $work_dir/${prefix}.hifiasm.log
# ## transform to fasta
# awk '/^S/{print ">"$2"\n"$3}' $work_dir/${prefix}.hifiasm.p_ctg.gfa > $work_dir/${prefix}.p_ctg.fa  ## chatgpt
awk '$1=="S" {printf ">%s\n%s\n", $2, $3} ' $work_dir/${prefix}.hifiasm.p_ctg.gfa > $work_dir/${prefix}.p_ctg.fa  ## rohan

## binning ggbin 

## find circular contigs  /shared/software/cobra-meta/latest/bin/cobra-meta
## scripts to parse the hifiasm output by rohan
/home/rohan/dev/pipeline/workflow/scripts/gg_rename_assembly.py \
  -i pbio-2857.29455.bc1003_BAK8A_OA--bc1003_BAK8A_OA.hifi_reads.bc1003_BAK8A_OA.ccs.trim.fq.gz.hifiasm_meta.p_ctg.gfa.fna \ 
  -o USGome_A.hifiasm.fna -s USGome_A

## identify MGE , virsorter2, virsorter1, vibrant, plasX,        unkown circular contigs possible MGE
# genomad end-to-end --relaxed --cleanup --enable-score-calibration \
# --threads $threads --sensitivity 7.0 --force-auto \
# $work_dir/${prefix}.p_ctg.fa \
# ${work_dir}_genomad/ \
# /groups/diamond/databases/genomad/v1.7/

# cat ${work_dir}_genomad/*/*plasmid_summary.tsv  ${work_dir}_genomad/*/*virus_summary.tsv > $work_dir/${prefix}.mge_summary.tsv


# ## map reads to assembly
# samtools faidx $work_dir/${prefix}.p_ctg.fa
# ~/smrtlink/pbmm2 align --preset CCS -j $threads $work_dir/${prefix}.p_ctg.fa $hifi_bam $work_dir/${prefix}.raw.bam
# samtools sort -T $work_dir/${prefix} -@ $threads -o $work_dir/${prefix}.align.bam $work_dir/${prefix}.raw.bam
# rm $work_dir/${prefix}.raw.bam
# samtools index $work_dir/${prefix}.align.bam
# /home/shuaiw//smrtlink/pbindex $work_dir/${prefix}.align.bam


/usr/bin/time -v -o $work_dir/${prefix}.methyl.time python /home/shuaiw/Methy/main.py \
  --work_dir $work_dir/${prefix}_methylation \
  --whole_bam $work_dir/${prefix}.align.bam \
  --whole_ref $work_dir/${prefix}.p_ctg.fa \
  --read_type hifi \
  --min_len 1000 \
  --max_NM 3000 \
  --min_cov 1 \
  --min_frac 0.4 \
  --min_score 30 \
  --min_sites 30 \
  --threads $threads \
  --plasmid_file $work_dir/${prefix}.mge_summary.tsv 


## based on bins and circular ones, clean up, checkM2 bins, separate the circular ones and the rest in the bin, 
## run each contig and each bin using checkM2, separate the MGE, waiting...

# gtdbtk classify_wf \
#   --genome_dir /home/shuaiw/methylation/data/borg/contigs/bins/ \
#   --out_dir /home/shuaiw/methylation/data/borg/contigs/GTDB \
#   --skip_ani_screen \
#   --cpus 64 \
#   -x fa

### prodigal-gv for predict faa, and gff, and prokka for rRNA, tRNA and KEGG, and dram for pfam, cazy,




