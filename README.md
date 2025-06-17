# Methy

## env requirement
```
conda activate /home/shuaiw/miniconda3/envs/methy3
```
## usage

```
python /home/shuaiw/Methy/main.py \
  --work_dir /home/shuaiw/methylation/data/borg/new_test12 \
  --whole_bam /home/shuaiw/methylation/data/borg/b_contigs/11.align.bam \
  --whole_ref /home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa \
  --read_type subreads 
```

## parameters
```
usage: main.py [-h] --whole_bam WHOLE_BAM --whole_ref WHOLE_REF --work_dir WORK_DIR [--maxAlignments MAXALIGNMENTS] [--read_type {subreads,hifi}] [--max_NM MAX_NM] [--min_len MIN_LEN] [--min_cov MIN_COV]
               [--kmer_mean_db KMER_MEAN_DB] [--kmer_num_db KMER_NUM_DB] [--no-clean] [--min_frac MIN_FRAC] [--min_sites MIN_SITES] [--min_score MIN_SCORE] [--plasmid_file PLASMID_FILE] [--threads THREADS]

Run methylation-based MGE-host linkage discovery pipeline.

options:
  -h, --help            show this help message and exit
  --whole_bam WHOLE_BAM
                        Input BAM file with kinetic data (HiFi or subreads).
  --whole_ref WHOLE_REF
                        Reference FASTA file for contigs.
  --work_dir WORK_DIR   Working directory for all output files.
  --maxAlignments MAXALIGNMENTS
                        Maximum number of alignments to process.
  --read_type {subreads,hifi}
                        Type of reads in BAM file.
  --max_NM MAX_NM       Maximum number of mismatches allowed (None to disable).
  --min_len MIN_LEN     Minimum contig length to process.
  --min_cov MIN_COV     Minimum read coverage required to retain a base.
  --kmer_mean_db KMER_MEAN_DB
                        Path to optional k-mer mean IPD database.
  --kmer_num_db KMER_NUM_DB
                        Path to optional k-mer count database.
  --no-clean            Disable cleaning step
  --min_frac MIN_FRAC   Minimum methylation fraction to retain a motif.
  --min_sites MIN_SITES
                        Minimum number of methylated sites per motif.
  --min_score MIN_SCORE
                        Minimum score for modification calling.
  --plasmid_file PLASMID_FILE
                        Optional plasmid FASTA file (set to 'NA' if not used).
  --threads THREADS     Number of threads to use for processing.
```




## Output interpretation

`figs/*.png`
| Subplot | Description |
| --- | --- |
|(0, 0)|Distribution of alignment coverage of all loci|
|(0, 1)|Distribution of raw IPD values|
|(1, 0)|Distribution of control IPD values|
|(1, 1)|Distribution of IPD ratios|

`gffs/*.reprocess.gff`
| Column | Description |
| --- | --- |
|1|contig name|
|2|method|
|3|methylation type|
|4|methylation position|
|5|methylation position|
|6|score|
|7|strand|
|8|.|
|9|additional annotation including coverage, local context, and IPDRtio of this locus|

`motifs/*motifs.csv`


`profiles/*.motifs.profile.csv`
| Column | Description |
| --- | --- |
|motifString| motif string|
|centerPos|Methylated locus on the string|
|for_loci_num|No. of motif sites on forward strand|
|for_modified_num|No. of methylated motif sites on forward strand|
|for_modified_ratio|Fraction of methylation in all motif sites on the forward strand|
|rev_loci_num|No. of motif sites on reverse strand|
|rev_modified_num|No. of methylated motif sites on reverse strand|
|rev_modified_ratio|Fraction of methylation in all motif sites on the reverse strand|
|motif_loci_num|No. of motif sites on both strand|
|motif_modified_num|No. of methylated motif sites on both strand|
|motif_modified_ratio|Fraction of methylation in all motif sites on the both strand|
|proportion|Fraction of methylation of this motif and and methylation of all motifs|

`hosts/`


`summary files`
| File | Description |
| --- | --- |
|motif_profile.csv| Methylation fraction of each motif on each contig|
|motif_cluster.csv|Contig clustering result based on motif profile|
|motif_heatmap.pdf| Visualization of motif profile|
|motif_cluster.pdf|Visualization of contig clusters|


## some notes
https://github.com/PacificBiosciences/pbcore/issues/127
```
  File "/home/shuaiw/miniconda3/envs/methy/lib/python3.12/site-packages/pbcore/io/align/BamIO.py", line 85, in _loadReadGroupInfo
    rgID = rgAsInt(rg["ID"])
           ^^^^^^^^^^^^^^^^^
  File "/home/shuaiw/miniconda3/envs/methy/lib/python3.12/site-packages/pbcore/io/align/_BamSupport.py", line 66, in rgAsInt
    return np.int32(int(rgIdString.split("/")[0], 16))
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
```
OverflowError: Python integer 3367457666 out of bounds for int32

```
pbcore
pip install git+https://github.com/PacificBiosciences/pbcore.git
numpy version

seaborn
scipy 
samtools
```
