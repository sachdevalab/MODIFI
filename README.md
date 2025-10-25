# 🧬 mGlu: DNA Modification Detection from PacBio SMRT Metagenomic Data

[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)
[![Python](https://img.shields.io/badge/Python-≥3.9-blue.svg)](#)
[![Platform](https://img.shields.io/badge/Platform-Linux-lightgrey.svg)](#)
[![PacBio](https://img.shields.io/badge/PacBio-SMRT-yellow.svg)](#)

**mGlu** is a software package for detecting DNA base modifications and inferring host–mobile genetic element (MGE) linkages from **PacBio SMRT metagenomic sequencing** data.  
It enables precise methylation calling, motif discovery, and host–MGE association in complex microbial communities.

---

## 🧩 Installation

### 1️⃣ Create and activate environment

```bash
conda env create -f env.yml
conda activate mGlu
```

### 2️⃣ Install dependencies

```bash
pip install git+https://github.com/PacificBiosciences/pbcore.git
```

> ⚠️ Note: Some PacBio BAMs may trigger an integer overflow error in `pbcore`.  
> See [PacificBiosciences/pbcore#127](https://github.com/PacificBiosciences/pbcore/issues/127) for details.

### 3️⃣ Compile C++ components

```bash
cd src && bash install.sh
```

### 4️⃣ Install SMRT Link Software

mGlu requires several tools from PacBio's SMRT Link software suite:

- **pbmotifmaker** - For motif discovery
- **pbmm2** - For sequence alignment (when using `--unaligned_bam`)
- **pbindex** - For BAM file indexing

**Installation options:**

1. **Install SMRT Link from PacBio** (smrtlink-release_25.3.0.273777+):
   Download and install from [PacBio SMRT Link](https://www.pacb.com/support/software-downloads/)

**Configuration:**
After installation, update the paths in `main.py` to point to your SMRT Link installation:
```python
motif_maker_bin = "/path/to/your/smrtlink/pbmotifmaker"
pbmm2_bin = "/path/to/your/smrtlink/pbmm2"
pbindex_bin = "/path/to/your/smrtlink/pbindex"
```

## 🧪 Test Installation

To verify that mGlu is working correctly, you can run the provided test dataset:

```bash
cd test/
sh test.sh
```

---

## ⚡ Quick Start

mGlu accepts both **aligned** and **unaligned** PacBio BAM files as input.

### Option A: Using pre-aligned BAM files

If you already have aligned BAM files:

```bash
python /home/shuaiw/Methy/main.py \
  --whole_bam /path/to/your/aligned.bam \
  --whole_ref /path/to/reference.fa \
  --work_dir /path/to/output \
  --read_type hifi
```

### Option B: Using unaligned BAM files (automatic alignment)

If you have unaligned BAM files, mGlu will automatically align them using pbmm2:

```bash
python /home/shuaiw/Methy/main.py \
  --unaligned_bam /path/to/your/unaligned.bam \
  --whole_ref /path/to/reference.fa \
  --work_dir /path/to/output \
  --read_type hifi
```

---



---

## ⚙️ Command-line Options

```bash
usage: main.py [-h] [--whole_bam WHOLE_BAM | --unaligned_bam UNALIGNED_BAM] 
               --whole_ref WHOLE_REF --work_dir WORK_DIR
               [--read_type {subreads,hifi}] [--max_NM MAX_NM] 
               [--min_len MIN_LEN] [--min_cov MIN_COV]
               [--kmer_mean_db KMER_MEAN_DB] [--kmer_num_db KMER_NUM_DB]
               [--min_frac MIN_FRAC] [--min_sites MIN_SITES]
               [--min_score MIN_SCORE] [--threads THREADS] [--clean]
```

### Input Options (choose one)

| Option | Description |
|--------|--------------|
| `--whole_bam` | Input **aligned** BAM file with kinetic data (HiFi or subreads). |
| `--unaligned_bam` | Input **unaligned** BAM file with kinetic data. Will be automatically aligned using pbmm2. |

### Required Parameters

| Option | Description |
|--------|--------------|
| `--whole_ref` | Reference FASTA file for contigs. |
| `--work_dir` | Working directory for all output files. |
| `--read_type` | Read type: `subreads` or `hifi`. Determines pbmm2 alignment preset (SUBREAD vs CCS). |

### Optional Parameters

| Option | Description | Default |
|--------|--------------|---------|
| `--max_NM` | Max number of mismatches allowed. | Auto (20M for HiFi, 10M for subreads) |
| `--min_len` | Minimum contig length to process. | 1000 |
| `--min_cov` | Minimum read coverage required per base. | 1 |
| `--min_ctg_cov` | Minimum read coverage required per contig. | 5 |
| `--kmer_mean_db` | Path to optional k-mer mean IPD database. | None |
| `--kmer_num_db` | Path to optional k-mer count database. | None |
| `--min_frac` | Minimum methylation fraction to retain a motif. | 0.4 |
| `--min_sites` | Minimum number of methylated sites per motif. | 30 |
| `--min_score` | Minimum score for modification calling. | 30 |
| `--threads` | Number of threads to use. | 64 |
| `--clean` | Remove intermediate files after completion. | False |

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
| `align_bam/aligned.bam` | Aligned BAM file (only created when using `--unaligned_bam`). |

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

