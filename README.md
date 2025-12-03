# 🧬 mGlu: DNA Modification Detection from PacBio SMRT Metagenomic Data

[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)
[![Python](https://img.shields.io/badge/Python-≥3.9-blue.svg)](#)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey.svg)](#)
[![PacBio](https://img.shields.io/badge/PacBio-SMRT-yellow.svg)](#)

**mGlu** is a software package for detecting DNA base modifications and inferring host–mobile genetic element (MGE) linkages from **PacBio SMRT metagenomic sequencing** data.  
It enables precise methylation calling, motif discovery, and host–MGE association in complex microbial communities.

---

## 🧩 Installation

### 1️⃣ Create and activate environment

```bash
conda create -n methy3 python=3.12
conda activate methy3
```

### 2️⃣ Install dependencies

```bash
pip install git+https://github.com/PacificBiosciences/pbcore.git
pip install numpy seaborn scipy
conda install -c bioconda samtools
```

> ⚠️ Note: Some PacBio BAMs may trigger an integer overflow error in `pbcore`.  
> See [PacificBiosciences/pbcore#127](https://github.com/PacificBiosciences/pbcore/issues/127) for details.

### 3️⃣ Add the `motifMaker` utility

```bash
add motifMaker
```
SMRT Link server software cannot be installed on Mac OS or Windows systems. So mGlu not support Mac OS and Windows.
---

## ⚡ Quick Start

### Step 1. Align PacBio reads to contigs

```bash
~/smrtlink/pbmm2 align --preset CCS -j $threads   $work_dir/${prefix}.p_ctg.fa $hifi_bam $work_dir/${prefix}.raw.bam

samtools sort -T $work_dir/${prefix} -@ $threads   -o $work_dir/${prefix}.align.bam $work_dir/${prefix}.raw.bam

rm $work_dir/${prefix}.raw.bam
samtools index $work_dir/${prefix}.align.bam
/home/shuaiw/smrtlink/pbindex $work_dir/${prefix}.align.bam
```

### Step 2. Run mGlu for methylation detection

```bash
python /home/shuaiw/Methy/main.py   --work_dir /home/shuaiw/methylation/data/borg/new_test12   --whole_bam /home/shuaiw/methylation/data/borg/b_contigs/11.align.bam   --whole_ref /home/shuaiw/methylation/data/borg/b_contigs/contigs/11.fa   --read_type subreads
```

---

## ⚙️ Command-line Options

```bash
usage: main.py [-h] --whole_bam WHOLE_BAM --whole_ref WHOLE_REF --work_dir WORK_DIR
               [--maxAlignments MAXALIGNMENTS] [--read_type {subreads,hifi}]
               [--max_NM MAX_NM] [--min_len MIN_LEN] [--min_cov MIN_COV]
               [--kmer_mean_db KMER_MEAN_DB] [--kmer_num_db KMER_NUM_DB]
               [--no-clean] [--min_frac MIN_FRAC] [--min_sites MIN_SITES]
               [--min_score MIN_SCORE] [--plasmid_file PLASMID_FILE] [--threads THREADS]
```

| Option | Description |
|--------|--------------|
| `--whole_bam` | Input BAM file with kinetic data (HiFi or subreads). |
| `--whole_ref` | Reference FASTA file for contigs. |
| `--work_dir` | Working directory for all output files. |
| `--read_type` | Read type: `subreads` or `hifi`. |
| `--max_NM` | Max number of mismatches allowed. |
| `--min_len` | Minimum contig length to process. |
| `--min_cov` | Minimum read coverage required per base. |
| `--kmer_mean_db` | Path to optional k-mer mean IPD database. |
| `--min_frac` | Minimum methylation fraction to retain a motif. |
| `--min_sites` | Minimum number of methylated sites per motif. |
| `--min_score` | Minimum score for modification calling. |
| `--threads` | Number of threads to use. |

---

## 📁 Output Overview

| File / Folder | Description |
|----------------|--------------|
| `motif_profile.csv` | Methylation fraction of each motif per contig. |
| `motif_heatmap.pdf` | Heatmap visualization of motif profiles. |
| `host_summary.csv` | Inferred best host for each MGE. |
| `gffs/*.reprocess.gff` | Base-level methylation annotations. |
| `motifs/*motifs.csv` | List of detected motifs and metrics. |
| `profiles/*.motifs.profile.csv` | Strand-specific methylation ratios. |
| `figs/*.png` | Quality-control and IPD distribution plots. |

---

### 📊 Example: GFF Annotation

| Column | Description |
|---------|-------------|
| 1 | Contig name |
| 2 | Method |
| 3 | Methylation type |
| 4–5 | Position |
| 6 | Score |
| 7 | Strand |
| 9 | Additional annotation (coverage, local context, IPD ratio) |

### 📈 Example: Figures (`figs/*.png`)

| Subplot | Description |
|----------|-------------|
| (0,0) | Alignment coverage distribution |
| (0,1) | Raw IPD value distribution |
| (1,0) | Control IPD distribution |
| (1,1) | IPD ratio distribution |

---

## 🧠 Control IPD (k-mer) Database

For isolated genomes or low-complexity metagenomes, provide an independent control dataset.  
The control IPD values are stored in:

```
control/control_db.up7.down3.mean.dat
control/control_db.up7.down3.num.dat
```

The same k-mer window size (up7/down3) must be used for both the sample and control database.

---

## 📊 MGE-Host Linkage Output Format

The MGE-host linkage prediction results contain the following columns:

| Column | Description |
|--------|-------------|
| `MGE` | Mobile genetic element (MGE) contig name |
| `MGE_len` | Length of the MGE contig (0 indicates circular element) |
| `host` | Predicted host contig name |
| `final_score` | Overall linkage confidence score (0-1, higher = more confident) |
| `specificity` | Specificity p-value for MGE-host association |
| `pvalue` | Statistical significance p-value for the linkage |
| `self_pvalue` | Self-comparison p-value (control) |
| `MGE_gc` | GC content of the MGE contig |
| `host_gc` | GC content of the host contig |
| `cos_sim` | Cosine similarity between MGE and host methylation profiles |
| `MGE_cov` | Sequencing coverage depth of the MGE |
| `host_cov` | Sequencing coverage depth of the host |
| `host_motif_num` | Number of methylation motifs detected in the host |
| `confidence` | Overall confidence score for the prediction |
| `motif_confidence` | Confidence score specific to motif-based evidence |
| `total_sites` | Total number of methylated sites used in the analysis |
| `motif_info` | Detailed motif information (format: motif:length:sites:coverage:score:position) |

---

## ⚠️ Notes

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

> *Please cite the corresponding publication once available.*

---

## 💬 Getting Help

For questions or bug reports:
- 📧 **Email:** [wshuai294@gmail.com](mailto:wshuai294@gmail.com)  
- 🐛 **GitHub Issues:** open an issue in this repository

We’ll respond as soon as possible.

