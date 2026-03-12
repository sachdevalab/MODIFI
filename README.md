# 🧬 MODIFI: DNA Modification Detection from PacBio SMRT Metagenomic Data

[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)
[![Python](https://img.shields.io/badge/Python-≥3.9-blue.svg)](#)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey.svg)](#)
[![PacBio](https://img.shields.io/badge/PacBio-SMRT-yellow.svg)](#)

**MODIFI** is a software package for detecting DNA base modifications and inferring host–mobile genetic element (MGE) linkages from **PacBio metagenomic sequencing** data. It enables precise modification calling, motif discovery, and host–MGE association in complex microbial communities, supporting both `subreads` and `HiFi` read types.

---

## 🧩 Installation

### 1️⃣ Clone repository and create environment

```bash
git clone https://github.com/wshuai294/MODIFI.git
cd MODIFI/
conda env create -n modifi -f env.yml
conda activate modifi
```

### 2️⃣ Install MODIFI via pip (recommended)

This will:
- Compile the C++ helper binary `get_control_IPD` (using `src/install.sh` under the hood)
- Install the full pipeline under the environment prefix
- Expose `modifi` / `MODIFI` commands on your `$PATH`

```bash
pip install .
```

### 3️⃣ Configure SMRT Link tools ⚠️ **OPTIONAL - Can be skipped**

> 💡 **You can skip this step!** MODIFI will automatically use tools from conda or fallback to built-in alternatives.

MODIFI requires three PacBio SMRT Link tools: `pbmotifmaker`, `pbmm2`, and `pbindex`.

**Configuration priority (in order):**

1. **Config file first** – If `config.yaml` exists, MODIFI will use the path specified there
2. **System PATH fallback** – If config.yaml is not found or incomplete, MODIFI checks system PATH
3. **MultiMotifMaker.jar/Conda** – Used as fallback for motif calling if `SMRT Link tools` is unavailable

**To configure via config.yaml (only if needed):**

Create or edit `config.yaml` in the MODIFI directory:

```yaml
smrtlink_bin: /path/to/smrtlink/private/bin/
```

> **Note:** If you have SMRT Link tools in your PATH (e.g., via conda), no configuration is needed. MODIFI will detect and use them automatically.

### 4️⃣ Verify installation

Test that the command-line entry point is available:

```bash
modifi --help
```

Or run the built-in test dataset:

```bash
cd test/
bash test.sh
```

### 5️⃣ [Optional] Setup for subreads

> ⚠️ **Note:** We recommend using HiFi reads instead of subreads when possible.

If you need to process subreads, create a separate environment with `pbcore` (requires Python 3.9 and numpy 1.22.4):

```bash
conda env create -n MODIFI_subreads -f subreads.env.yml
conda activate MODIFI_subreads
pip install git+https://github.com/PacificBiosciences/pbcore.git
```

---

## ⚡ Quick Start

### Input Requirements

MODIFI supports both `HiFi` and `subreads`. Your PacBio BAM file **must contain kinetics (IPD) tags**:

| Read Type | Required IPD Tags | Description | 
|-----------|-------------------|-------------| 
| HiFi | `fi`, `ri` | Forward/Reverse IPD (codec V1) |
| Subreads | `ip` | IPD (raw frames or codec V1) |

> ⚠️ **Important:** 
> - FASTQ files do not contain kinetics information and **cannot be used**
> - The BAM file stores kinetic information; it does not need to be pre-aligned
> - See [PacBio BAM format specifications](https://pacbiofileformats.readthedocs.io/en/11.0/BAM.html) for details

### Basic Usage

Run MODIFI with unaligned BAM containing kinetics:

```bash
modifi \
    --work_dir /path/to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref assembly.fasta \
    --read_type hifi
```

For all options:
```bash
modifi --help
```

---

## 📖 Detailed Usage

### Input BAM Options

> 💡 **Tip:** Once motif detection has finished (i.e. `profiles/` exists in `--work_dir`), you can re-run **linkage only** using the `modifi-linkage` helper instead of the full pipeline (see [Linkage-only mode](#linkage-only-mode)).

#### Option 1: Unaligned BAM (recommended)

MODIFI will automatically align reads using pbmm2:

```bash
modifi \
    --work_dir /path/to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref assembly.fasta \
    --read_type hifi
```

#### Option 2: Pre-aligned BAM

If already aligned with pbmm2 (with kinetics preserved):

```bash
modifi \
    --work_dir /path/to/output \
    --whole_bam aligned.sorted.bam \
    --whole_ref assembly.fasta \
    --read_type hifi
```

#### Option 3: Subreads

For subread data (requires `MODIFI_subreads` environment):

```bash
modifi \
    --work_dir /path/to/output \
    --unaligned_bam raw.subreads.bam \
    --whole_ref assembly.fasta \
    --read_type subreads
```

### MGE-Host Linkage Inference

Provide an MGE table to automatically infer MGE-host linkages:

```bash
modifi \
    --work_dir /path/to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref assembly.fasta \
    --read_type hifi \
    --mge_file MGE_list.tab
```

**MGE table format:** Tab-separated file with at least a `seq_name` column. Can be output from [geNomad](https://github.com/apcamargo/genomad) or manually created:

```tsv
seq_name	length	topology	n_genes	genetic_code	plasmid_score	fdr
ERR6535514_2_L	6291	DTR	8	11	0.9997	0.0003
...
```

### Incorporating Binning Results

Assign MGEs to host bins using binning results:

```bash
modifi \
    --work_dir /path/to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref assembly.fasta \
    --read_type hifi \
    --mge_file MGE_list.tab \
    --bin_file binning.tab
```

**Bin file format:** Tab-separated, two columns (no header): contig name and bin name

```tsv
scaffold_53	maxbin.017
scaffold_66	maxbin.017
scaffold_73	maxbin.017
...
```

### Using Control Databases

For isolate genomes or low-complexity metagenomes, use a control database for better normalization.

#### Pre-built database:

```bash
modifi \
    --work_dir /path/to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref assembly.fasta \
    --read_type hifi \
    --kmer_mean_db control/control_db.up7.down3.mean.dat \
    --kmer_num_db control/control_db.up7.down3.num.dat
```

#### Custom database:

Use control database from your own high-complexity metagenome runs (automatically generated in `<work_dir>/control/`):

```
control/control_db.up7.down3.mean.dat
control/control_db.up7.down3.num.dat
```


---

---

## 🔬 Complete Metagenomics Workflow

End-to-end example from MGE prediction to host linkage inference.

**Prerequisites:** 
- Raw PacBio BAM with kinetics
- Metagenomic assembly (FASTA)

### Step 1: Predict MGEs with geNomad

```bash
genomad end-to-end \
    --relaxed \
    --cleanup \
    --enable-score-calibration \
    --threads 32 \
    --sensitivity 7.0 \
    --force-auto \
    assembly.fasta \
    output/genomad/ \
    /path/to/genomad_db
```

### Step 2: Combine virus and plasmid predictions

```bash
cat output/genomad/assembly_summary/assembly_virus_summary.tsv \
    output/genomad/assembly_summary/assembly_plasmid_summary.tsv \
    > all_MGEs.tsv
```

### Step 3: Run MODIFI

```bash
modifi \
    --work_dir output/MODIFI \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref assembly.fasta \
    --read_type hifi \
    --mge_file all_MGEs.tsv \
    --threads 32
```

> **Note:** Proviruses are automatically excluded from analysis.

---

### Linkage-only mode

If you have already run MODIFI and generated motif profiles in a given `--work_dir`, you can re-run just the **MGE-host linkage** step (for example, with a different MGE table or binning file) without re-running modification detection:

```bash
modifi-linkage \
    --work_dir /path/to/output \
    --whole_ref assembly.fasta \
    --mge_file all_MGEs.tsv \
    --bin_file binning.tab \
    --threads 32
```

This command:
- Reads existing motif profiles from `<work_dir>/profiles/`
- Writes host linkage results into `<work_dir>/hosts/`
- Uses the same linkage logic as the `host` step in the full pipeline.

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

## � Output Files

### Main Output Files

| File | Description |
|------|-------------|
| `motif_profile.csv` | Methylation fraction matrix: motifs × contigs |
| `host_summary.csv` | Best predicted host for each MGE |
| `mean_depth.csv` | Sequencing depth and length per contig |
| `summary.csv` | Overall modification detection statistics |
| `motif_heatmap.pdf` | Heatmap of motif profiles (empty for large datasets) |

### Per-Genome Results

| Folder | Contents |
|--------|----------|
| `gffs/` | `*.reprocess.gff` – Base-level methylation annotations |
| `motifs/` | `*.motifs.csv` – Detected motifs and metrics |
| `profiles/` | `*.motifs.profile.csv` – Strand-specific methylation ratios |
| `hosts/` | `*.host_prediction.csv`  – Linkage scores for all candidate hosts |
| `hosts/` | `*.motif_data.csv` – Detailed motif data for scoring |

### Optional Output

| Folder | Contents | Enabled by |
|--------|----------|------------|
| `RM_systems/` | RM gene annotations | `--annotate_rm` |
| `figs/` | IPD distribution plots | `--visu_ipd` |
| `ipd_ratio/` | IPD ratio line plots | `--detect_misassembly` |

---

### Output Format: Motif File

**File:** `motifs/*.motifs.csv`

| Column | Description |
|--------|-------------|
| `motifString` | DNA motif sequence (e.g., GATC, TCGAG) |
| `centerPos` | Modified base position within motif (1-indexed) |
| `modificationType` | Modification type (e.g., modified_base) |
| `fraction` | Fraction of modified sites (0–1) |
| `nDetected` | Number of modified sites detected |
| `nGenome` | Total motif sites in genome |
| `groupTag` | Motif group ID for clustering |
| `partnerMotifString` | Complementary motif (if paired) |
| `meanScore` | Mean modification score |
| `meanIpdRatio` | Mean IPD ratio |
| `meanCoverage` | Mean coverage at motif sites |
| `objectiveScore` | Overall detection confidence score |


### Output Format: MGE-Host Linkage

**File:** `host_summary.csv`

| Column | Description |
|--------|-------------|
| `MGE` | MGE contig name |
| `MGE_len` | MGE length (bp) |
| `host` | Predicted host contig/bin |
| `final_score` | Linkage confidence (0–1, higher is better) |
| `specificity` | Association specificity |
| `MGE_gc` | MGE GC content (%) |
| `host_gc` | Host GC content (%) |
| `cos_sim` | 4-mer frequency cosine similarity |
| `MGE_cov` | MGE sequencing depth |
| `host_cov` | Host sequencing depth |
| `host_motif_num` | Number of motifs in host |
| `total_sites` | Total modified sites analyzed |
| `motif_info` | Detailed motif data (format: motif:length:sites:cov:score:pos) |

> **Filtering recommendations:** 
> - High confidence: `final_score > 0.5` AND `specificity < 0.01`
> - Medium confidence: `final_score > 0.3` AND `specificity < 0.05`

### Output Format: GFF Annotation

**Files:** `gffs/*.reprocess.gff`

Standard GFF3 format with 9 columns:

| Column | Content |
|--------|----------|
| 1 | Contig name |
| 2 | Method (MODIFI) |
| 3 | Feature type (modified_base) |
| 4–5 | Position (start–end, 1-indexed) |
| 6 | Modification score |
| 7 | Strand (+/-) |
| 8 | Phase (not used) |
| 9 | Attributes (coverage, context, IPD ratio) |

### Quality Control Plots

**Files:** `figs/*.png` (generated with `--visu_ipd`)

2×2 subplot grid showing IPD normalization:

| Position | Content |
|----------|----------|
| (0,0) | Alignment coverage distribution |
| (0,1) | Raw IPD value distribution |
| (1,0) | Control IPD distribution |
| (1,1) | Normalized IPD ratio distribution |

---

## 🧠 Control Database

For **isolate genomes** or **low-complexity metagenomes**, using a control database improves IPD normalization and modification calling accuracy.

### Pre-built Database

MODIFI includes a pre-built control database:

```
control/control_db.up7.down3.mean.dat
control/control_db.up7.down3.num.dat
```

Use with `--kmer_mean_db` and `--kmer_num_db` flags (see [Using Control Databases](#using-control-databases)).

### Custom Database

For better results, generate a control database from your own **high-complexity metagenomes**:

1. Run MODIFI on high-complexity metagenomic samples
2. Control files are automatically generated in `<work_dir>/control/`
3. Reuse these files for isolate/low-complexity samples

### Important Notes

- The k-mer window size (e.g., `up7.down3`) **must match** between sample and control
- Control database accumulates more accurate statistics with larger/more diverse datasets
- Recommended for: isolates, enriched cultures, low-diversity communities

---

## ⚠️ Important Notes

### Performance
- Coverage **>500×** is automatically subsampled for computational efficiency
- For large metagenomes, use `--threads` to parallelize processing
- `--segment` flag increases recall for low-depth contigs but requires more time


### Best Practices
- Use **HiFi reads** instead of subreads when possible (better accuracy, faster)
- Use control database for isolate genomes and low-complexity samples
- Filter MGE-host linkages: `final_score > 0.5` and `specificity < 0.01`

---

## 📚 Citation

If you use MODIFI in your research, please cite:

> *Manuscript in preparation.*

---

## 💬 Getting Help

For questions or bug reports:
- 📧 **Email:** [wshuai294@gmail.com](mailto:wshuai294@gmail.com)  
- 🐛 **GitHub Issues:** open an issue in this repository

We’ll respond as soon as possible.

