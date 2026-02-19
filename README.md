# 🧬 mGlu: DNA Modification Detection from PacBio SMRT Metagenomic Data

[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)
[![Python](https://img.shields.io/badge/Python-≥3.9-blue.svg)](#)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey.svg)](#)
[![PacBio](https://img.shields.io/badge/PacBio-SMRT-yellow.svg)](#)

**mGlu** is a software package for detecting DNA base modifications and inferring host–mobile genetic element (MGE) linkages from **PacBio metagenomic sequencing** data. It enables precise modification calling, motif discovery, and host–MGE association in complex microbial communities, supporting both `subreads` and `HiFi` read types.

---

## 🧩 Installation

### 1️⃣ Create and activate environment

```bash
git clone https://github.com/wshuai294/mGlu.git
cd mGlu/
conda env create -n mglu -f env.yml
conda activate mglu
```

### 2️⃣ Add the `SMRT Link` path
mGlu requires three tools from SMRT Link: `pbmotifmaker`, `pbmm2`, and `pbindex`. If they are in the system path, mGlu will automatically detect them. If not, you can set the path for SMRT Link tools in `mGlu/config.yaml` as follows:

```yaml
smrtlink_bin: /path/to/your/smrtlink/
```

### 3️⃣ Test installation
Verify that mGlu is installed correctly by running the test suite:

```bash
cd test/
bash test.sh
```

### 4️⃣ Additional module required for subreads
Install the `pbcore` module only if you want to handle `subreads`.

```bash
pip install git+https://github.com/PacificBiosciences/pbcore.git
```



---

## ⚡ Quick Start

### Step 1: Input read type requirements
mGlu supports both `subreads` and `HiFi` read types. Ensure your input PacBio BAM file contains IPD tags. This should be a **BAM** file used to store kinetics information, not necessarily an alignment BAM. The required IPD tags are:

| Read Type | IPD Tag (kinetics) | Description | 
|-----------|--------------------|--------------| 
| hifi | `fi` | Forward IPD (codec V1) |
| hifi | `ri` | Reverse IPD (codec V1) |
| subreads | `ip` | IPD (raw frames or codec V1) |

> ⚠️ **Note:** FASTQ files do not contain kinetics information and cannot be used.

For more information about PacBio BAM format specifications, see: https://pacbiofileformats.readthedocs.io/en/11.0/BAM.html

### Step 2: Run mGlu for modification detection
Example command:

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type hifi 
```

For more information, run:
```bash
python mGlu/main.py --help
```

--- 
## 📖 Detailed Usage Guidelines

### Using unaligned BAM files
If your BAM file contains unaligned reads, mGlu will automatically align them using pbmm2. Just run

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type hifi 
```

### Using pre-aligned BAM files
If you have already aligned the reads using pbmm2, use the `--whole_bam` option:

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --whole_bam pbmm2.aligned.sorted.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type hifi 
```

### Using subreads
For subread data, specify `--read_type subreads`:

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --unaligned_bam raw.subreads.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type subreads 
```

### MGE-host linkage inference
If you have an MGE table file, provide it using `--mge_file`, and mGlu will automatically infer MGE-host linkages:

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type hifi \
    --mge_file MGE_list.tab
```

The MGE table file can be the output from [geNomad](https://github.com/apcamargo/genomad) or generated manually. It should be a tab-separated file with at least one column with the header `seq_name`. Here is an example:

```
seq_name	length	topology	n_genes	genetic_code	plasmid_score	fdr	n_hallmarks	marker_enrichment	conjugation_genes	amr_genes
ERR6535514_2_L	6291	DTR	8	11	0.9997	0.0003	2	1.8745	NA	NA
...
```

### Assigning MGEs to host bins
If you have binning results and want to assign MGEs to host bins, provide the bin file using `--bin_file`:

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type hifi \
    --mge_file MGE_list.tab \
    --bin_file my_bin.tab
```

The `bin_file` should be a tab-separated file with two columns (no header required): contig name and bin name. Example:

```
RuReacBro_20230708_10_0h_50ppm_r1_scaffold_53   RuReacBro_20230708_10_0h_50ppm_r1_maxbin.017_rmcirc
RuReacBro_20230708_10_0h_50ppm_r1_scaffold_66   RuReacBro_20230708_10_0h_50ppm_r1_maxbin.017_rmcirc
RuReacBro_20230708_10_0h_50ppm_r1_scaffold_73   RuReacBro_20230708_10_0h_50ppm_r1_maxbin.017_rmcirc
RuReacBro_20230708_10_0h_50ppm_r1_scaffold_77   RuReacBro_20230708_10_0h_50ppm_r1_maxbin.017_rmcirc
RuReacBro_20230708_10_0h_50ppm_r1_scaffold_84   RuReacBro_20230708_10_0h_50ppm_r1_maxbin.017_rmcirc
...
```

### Using control databases
For isolate genomes or low-complexity metagenomes, use our pre-built control database:

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type hifi \
    --kmer_mean_db mGlu/control/control_db.up7.down3.mean.dat \
    --kmer_num_db mGlu/control/control_db.up7.down3.num.dat
```

Alternatively, you can use the control database generated from your own high-complexity metagenomes. These files will be stored at:

```
high_complexity_mGlu_output/control/control_db.up7.down3.mean.dat 
high_complexity_mGlu_output/control/control_db.up7.down3.num.dat 
```


---

## 🔬 Standard Metagenomics Analysis Workflow

This example demonstrates a complete workflow from MGE prediction to host linkage inference.

**Prerequisites:** Raw PacBio BAM file and metagenomic assembly.

### Step 1: Predict MGEs using geNomad

```bash
genomad end-to-end --relaxed --cleanup --enable-score-calibration \
    --threads <threads> --sensitivity 7.0 --force-auto \
    <assembly.fasta> \
    <output_dir>/Genomad/ \
    <genomad_database_path>
```

### Step 2: Merge geNomad virus and plasmid results

```bash
cat <output_dir>/Genomad/sample_summary/sample_virus_summary.tsv \
    <output_dir>/Genomad/sample_summary/sample_plasmid_summary.tsv \
    > all_MGEs.tsv
```

### Step 3: Run mGlu for modification detection and MGE-host linkage inference

```bash
python mGlu/main.py \
    --work_dir /path-to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref metagenomic_assembly.fasta \
    --read_type hifi \
    --mge_file all_MGEs.tsv
```

---

## ⚙️ Command-line Options

```bash
usage: main.py [-h] [-v] [--whole_bam WHOLE_BAM] [--unaligned_bam UNALIGNED_BAM] --whole_ref WHOLE_REF --work_dir WORK_DIR
               [--read_type {subreads,hifi}] [--min_iden MIN_IDEN] [--min_len MIN_LEN] [--min_cov MIN_COV] [--min_ctg_cov MIN_CTG_COV]
               [--kmer_mean_db KMER_MEAN_DB] [--kmer_num_db KMER_NUM_DB] [--clean] [--segment] [--min_frac MIN_FRAC]
               [--min_sites MIN_SITES] [--min_score MIN_SCORE] [--mge_file MGE_FILE] [--bin_file BIN_FILE] [--threads THREADS] [--up UP]
               [--down DOWN] [--detect_misassembly] [--visu_ipd] [--binning] [--annotate_rm] [--rm_gene_file RM_GENE_FILE]
               [--run_steps {split,load,control,compare,motif,profile,merge,host,anno} [{split,load,control,compare,motif,profile,merge,host,anno} ...]]

Run methylation-based MGE-host linkage discovery pipeline.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --whole_bam WHOLE_BAM
                        Input aligned BAM file with kinetic data (HiFi or subreads). Use this for pre-aligned BAM files. Must be aligned
                        by pbmm2 to ensure kinetic tags are present. (default: None)
  --unaligned_bam UNALIGNED_BAM
                        Input unaligned BAM file with kinetic data (HiFi or subreads). Will be aligned using pbmm2. (default: None)
  --whole_ref WHOLE_REF
                        Reference FASTA file for contigs. (default: None)
  --work_dir WORK_DIR   Working directory for all output files. (default: None)
  --read_type {subreads,hifi}
                        Type of reads in BAM file. (default: hifi)
  --min_iden MIN_IDEN   Minimum identity allowed for read alignment (default: 0.97 for hifi, 0.85 for subreads). (default: None)
  --min_len MIN_LEN     Minimum contig length to process. (default: 1000)
  --min_cov MIN_COV     Minimum read coverage required to retain a base. (default: 1)
  --min_ctg_cov MIN_CTG_COV
                        Minimum read coverage required to retain a contig. (default: 5)
  --kmer_mean_db KMER_MEAN_DB
                        Path to optional k-mer mean IPD database. (default: None)
  --kmer_num_db KMER_NUM_DB
                        Path to optional k-mer count database. (default: None)
  --clean               Enable cleaning step (default: False)
  --segment             Enable segmentation of the contigs by depth to increase recall for low-depth contigs (costs more time).
                        (default: False)
  --min_frac MIN_FRAC   Minimum methylation fraction to retain a motif. (default: 0.4)
  --min_sites MIN_SITES
                        Minimum number of methylated sites per motif. (default: 30)
  --min_score MIN_SCORE
                        Minimum score for modification calling. (default: 30)
  --mge_file MGE_FILE   MGE table file (sep by tab), can be output of geNomad, with at least one column with header: seq_name. (default:
                        NA)
  --bin_file BIN_FILE   Path to the binning file containing contig-to-bin mappings. (default: None)
  --threads THREADS     Number of threads to use for processing. (default: 64)
  --up UP               Number of upstream bases to consider for k-mer analysis. (default: 7)
  --down DOWN           Number of downstream bases to consider for k-mer analysis. (default: 3)
  --detect_misassembly  Enable detection of misassembly in the pipeline. (default: False)
  --visu_ipd            Enable visualization of IPD distribution. (default: False)
  --binning             Enable binning based on methylation (in testing). (default: False)
  --annotate_rm         Enable RM system annotation, MicrobeMod should be installed. (default: False)
  --rm_gene_file RM_GENE_FILE
                        RM gene annotation file by MicrobeMod, with suffix .rm.genes.tsv (only for testing) (default: None)
  --run_steps {split,load,control,compare,motif,profile,merge,host,anno} [{split,load,control,compare,motif,profile,merge,host,anno} ...]
                        Steps to run in the pipeline (default: all), sep by space. Only for easy testing. (default: ['split', 'load',
                        'control', 'compare', 'motif', 'profile', 'merge', 'host'])
```



---

## 📁 Output Overview

| File / Folder | Description |
|----------------|--------------|
| `motif_profile.csv` | Methylation fraction of each motif per contig. |
| `motif_heatmap.pdf` | Heatmap visualization of motif profiles (may be empty for large datasets). |
| `host_summary.csv` | Inferred best host for each MGE. |
| `gffs/*.reprocess.gff` | Base-level methylation annotations for each genome. |
| `motifs/*.motifs.csv` | List of detected motifs and metrics for each genome. |
| `profiles/*.motifs.profile.csv` | Strand-specific methylation ratios for each genome. |
| `hosts/*.host_prediction.csv` | Linkage scores for all candidate hosts of each MGE. |
| `hosts/*.motif_data.csv` | Detailed motif data for linkage score calculation of each MGE. |
| `mean_depth.csv` | Sequencing depth and length for each contig. |
| `summary.csv` | Summary statistics of modification detection across all genomes. |
| `RM_systems/*` | RM (Restriction-Modification) gene annotation results. |
| `figs/*.png` | Quality-control and IPD distribution plots (generated with `--visu_ipd`). |
| `ipd_ratio/*.pdf` | Line plots of IPD ratio at each motif site along the genome (generated with `--detect_misassembly`). |

---

### 📊 Motif File Format (`motifs/*.motifs.csv`)

| Column | Description |
|--------|-------------|
| `motifString` | DNA sequence motif pattern (e.g., TCGAG) |
| `centerPos` | Position of the modified base within the motif (1-indexed) |
| `modificationType` | Type of modification detected (e.g., modified_base) |
| `fraction` | Fraction of motif sites that are modified (0-1) |
| `nDetected` | Number of modified sites detected for this motif |
| `nGenome` | Total number of motif sites present in the genome |
| `groupTag` | Motif group identifier for clustering related motifs |
| `partnerMotifString` | Complementary or related motif on opposite strand (if applicable) |
| `meanScore` | Mean modification score across all detected sites |
| `meanIpdRatio` | Mean IPD ratio across all detected sites |
| `meanCoverage` | Mean sequencing coverage across all motif sites |
| `objectiveScore` | Overall confidence score for motif detection |


### 📊 MGE-Host Linkage Output Format

The MGE-host linkage prediction results contain the following columns:

| Column | Description |
|--------|-------------|
| `MGE` | Mobile genetic element (MGE) contig name |
| `MGE_len` | Length of the MGE contig |
| `host` | Predicted host contig (bin) name |
| `final_score` | Overall linkage confidence score (0-1, higher = more confident) |
| `specificity` | Specificity for MGE-host association |
| `pvalue` | Statistical significance p-value for the linkage |
| `self_pvalue` | Self-comparison p-value (control) |
| `MGE_gc` | GC content of the MGE contig |
| `host_gc` | GC content of the host contig |
| `cos_sim` | Cosine similarity between MGE and host 4-mer frequency |
| `MGE_cov` | Sequencing coverage depth of the MGE |
| `host_cov` | Sequencing coverage depth of the host |
| `host_motif_num` | Number of modification motifs detected in the host |
| `confidence` | Overall confidence score for the prediction |
| `motif_confidence` | Confidence score specific to motif-based evidence |
| `total_sites` | Total number of modified sites used in the analysis |
| `motif_info` | Detailed motif information (format: motif:length:sites:coverage:score:position) |

> **Recommended filtering criteria:** For confident MGE-host linkages, use `final_score > 0.5` and `specificity < 0.01`.

### 📊 DNA Modification GFF Annotation

| Column | Description |
|---------|-------------|
| 1 | Contig name |
| 2 | Method |
| 3 | Modification type (placeholder, not used) |
| 4–5 | Position (start–end) |
| 6 | Score |
| 7 | Strand |
| 8 | Frame (not used) |
| 9 | Additional annotation (coverage, local context, IPD ratio) |

### 📈 Illustration of IPD normalization (`figs/*.png`)

| Subplot | Description |
|----------|-------------|
| (0,0) | Alignment coverage distribution |
| (0,1) | Raw IPD value distribution |
| (1,0) | Control IPD distribution |
| (1,1) | IPD ratio distribution |

---

## 🧠 Control IPD (k-mer) Database

For isolated genomes or low-complexity metagenomes, an independent control dataset is recommended.  
The control IPD values are stored in:

```
control/control_db.up7.down3.mean.dat
control/control_db.up7.down3.num.dat
```

> **Note:** The k-mer window size (e.g., up7/down3) must match between your sample and the control database.




---

## ⚠️ Important Notes

- Coverage > **500×** is automatically subsampled for computational efficiency.  
- If you encounter:
  ```text
  OverflowError: Python integer 3367457666 out of bounds for int32
  ```
  reinstall pbcore via:
  ```bash
  pip install --force-reinstall git+https://github.com/PacificBiosciences/pbcore.git
  ```

---

## 📚 Citation

If you use mGlu in your research, please cite:

> *Manuscript in preparation.*

---

## 💬 Getting Help

For questions or bug reports:
- 📧 **Email:** [wshuai294@gmail.com](mailto:wshuai294@gmail.com)  
- 🐛 **GitHub Issues:** open an issue in this repository

We’ll respond as soon as possible.

