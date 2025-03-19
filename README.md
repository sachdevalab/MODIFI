# Methy

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
|9|additional annotation including coverage, local contex, and IPDRtio of this locus|

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

`summary files`
| File | Description |
| --- | --- |
|motif_profile.csv| Methylation fraction of each motif on each contig|
|motif_cluster.csv|Contig clustering result based on motif profile|
|motif_heatmap.pdf| Visualization of motif profile|
|motif_cluster.pdf|Visualization of contig clusters|
