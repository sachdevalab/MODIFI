# 🧬 MODIFI: DNA Modification Detection from PacBio SMRT Metagenomic Data

[![License](https://img.shields.io/badge/License-MIT-green.svg)](./LICENSE)
[![Python](https://img.shields.io/badge/Python-≥3.9-blue.svg)](#)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey.svg)](#)
[![PacBio](https://img.shields.io/badge/PacBio-SMRT-yellow.svg)](#)

**MODIFI** is a software package for detecting DNA base modifications and inferring host–mobile genetic element (MGE) linkages from **PacBio metagenomic sequencing** data. It enables precise modification calling, motif discovery, and host–MGE association in complex microbial communities, supporting both `subreads` and `HiFi` read types.

### Table of contents

- [How it works](#how-it-works)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Detailed Usage](#-detailed-usage)
- [Complete Metagenomics Workflow](#-complete-metagenomics-workflow)
- [Linkage-only mode](#linkage-only-mode)
- [Command-line options](#-command-line-options)
- [Output Files](#-output-files)
- [Control Database](#-control-database)
- [Important Notes](#-important-notes)
- [Citation](#-citation)
- [Getting Help](#-getting-help)

---

## How it works

MODIFI expects a **PacBio BAM with IPD kinetics** and a **reference FASTA** (metagenome assembly or isolate). If you pass an unaligned BAM, reads are aligned with pbmm2 while preserving kinetic tags. The pipeline splits work by contig, optionally accumulates **control k-mer IPD** statistics (from the same run or external databases), and **normalizes** sample IPDs against those controls. It then **calls modifications**, **discovers motifs**, and builds **per-contig methylation profiles**. With an MGE table (e.g. from geNomad) and optionally a **binning** file, it **merges** motif evidence and scores **host–MGE linkages**. Typical outputs are motif calls, strand-specific profiles, base-level GFF annotations, and `host_summary.csv` when linkage is enabled.

---

## 🧩 Installation

### 1️⃣ Clone repository and create environment

```bash
git clone https://github.com/wshuai294/MODIFI.git
cd MODIFI/
conda env create -n modifi -f env.yml
conda activate modifi
```

### 2️⃣ Install MODIFI via pip 

This will:
- Compile the C++ helper binary `get_control_IPD` (using `src/install.sh` under the hood)
- Install the full pipeline under the environment prefix
- Expose `modifi` / `MODIFI` commands on your `$PATH`

```bash
pip install .
```

### 3️⃣ Configure SMRT Link tools **OPTIONAL - Can be skipped**

> 💡 **You can skip this step!** MODIFI will automatically use tools from conda or fallback to built-in alternatives.

MODIFI requires three PacBio SMRT Link tools: `pbmotifmaker`, `pbmm2`, and `pbindex`.

**Configuration priority (in order):**

1. **Config file first** – If `smrt.config.yaml` exists, MODIFI will use the path specified there
2. **System PATH fallback** – If config.yaml is not found or incomplete, MODIFI checks system PATH
3. **MultiMotifMaker.jar/Conda** – Used as fallback if `SMRT Link tools` is unavailable

**To configure via config.yaml (only if needed):**

Create or edit `smrt.config.yaml` in the MODIFI directory:

```yaml
smrtlink_bin: /path/to/smrtlink/private/bin/
```

> **Note:** If you have SMRT Link tools in your PATH (e.g., via conda), no configuration is needed. MODIFI will detect and use them automatically.

### 4️⃣ Verify installation

Test that the command-line entry point is available:

```bash
modifi --help
```

And run the built-in test dataset:

```bash
cd test/hifi/
bash test_hifi.sh
```

### 5️⃣ [Optional] Setup for subreads

If you need to process subreads, create a separate environment with `pbcore` (requires Python 3.9 and numpy 1.22.4):

```bash
conda env create -n MODIFI_subreads -f subreads.env.yml
conda activate MODIFI_subreads
pip install git+https://github.com/PacificBiosciences/pbcore.git
pip install .
```

---

## ⚡ Quick Start

### Input Requirements

MODIFI supports both `HiFi reads` and `subreads`. Your PacBio BAM file **must contain kinetics (IPD) tags**:

| Read Type | Required IPD Tags | Description | 
|-----------|-------------------|-------------| 
| HiFi | `fi`, `ri` | Forward/Reverse IPD (codec V1) |
| Subreads | `ip` | IPD (raw frames or codec V1) |

> ⚠️ **Important:** 
> - FASTQ files do not contain kinetics information and **cannot be used**
> - The BAM file stores kinetic information; it does not need to be pre-aligned
> - See [PacBio BAM format specifications](https://pacbiofileformats.readthedocs.io/en/11.0/BAM.html) for details

### Example

Run MODIFI with unaligned BAM containing kinetics:

```bash
modifi \
    --work_dir /path/to/output \
    --unaligned_bam raw.hifi_reads.bam \
    --whole_ref assembly.fasta \
    --read_type hifi
```

---

## 📖 Detailed Usage

### Input BAM Options

#### Option 1: Unaligned BAM 

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
> 💡 **Tip:** Once motif detection has finished, you can re-run **linkage only** using the `modifi-linkage` helper instead of the full pipeline (see [Linkage-only mode](#linkage-only-mode)).

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

---

## ⚙️ Command-line options

For the full, up-to-date list of arguments and defaults:

```bash
modifi --help
modifi-linkage --help
```

Commonly tuned options (`modifi`):

| Option | Role |
|--------|------|
| `--work_dir` | Output directory for all pipeline artifacts |
| `--unaligned_bam` / `--whole_bam` | Unaligned (aligned with pbmm2 in-pipeline) or pre-aligned BAM with kinetics |
| `--whole_ref` | Reference FASTA (contigs) |
| `--read_type` | `hifi` (default) or `subreads` |
| `--threads` | Parallelism (default 64) |
| `--min_iden` | Minimum alignment identity (default 0.97 for HiFi, 0.85 for subreads if unset) |
| `--mge_file` / `--bin_file` | MGE table (e.g. geNomad) and optional contig→bin map for linkage |
| `--kmer_mean_db` / `--kmer_num_db` | Optional control k-mer IPD databases; window must match `--up` / `--down` |
| `--up` / `--down` | Upstream / downstream k-mer flank length around the modified base (defaults 7 / 3) |
| `--min_len`, `--min_cov`, `--min_ctg_cov` | Contig and site coverage filters |
| `--min_frac`, `--min_sites`, `--min_score` | Motif retention and modification calling thresholds |
| `--segment` | Depth-based contig segmentation (more recall, more runtime) |
| `--run_steps` | Restrict to a subset of steps: `split`, `load`, `control`, `compare`, `motif`, `profile`, `merge`, `host`, `anno` |
| `--clean` | Enable cleaning step |
| `--visu_ipd` | Write IPD distribution figures under `figs/` |
| `--detect_misassembly` | Misassembly-related IPD ratio outputs |
| `--annotate_rm` / `--rm_gene_file` | RM-system annotation (MicrobeMod; testing) |
| `--binning` | Methylation-based binning (testing) |


---

## 📁 Output Files

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

### Resource requirements (examples for reference)
In 59 metagenomics from nine habitats, the mean wall-clock time was 4.5 hours (range: 0.2–19 hours, std: 3.9), while the mean CPU time was 120 hours (range: 1–1,161 hours, std: 215), with a mean peak memory usage of 18 GB (range: 0.6–62 GB, std: 14). 

### Performance
- Coverage **>500×** is automatically subsampled for computational efficiency
- For large metagenomes, use `--threads` to parallelize processing
- `--segment` flag increases recall for low-depth contigs but requires more time

### Benchmarks and reproducibility

Driver scripts and plotting code used for analyses in development live under [`benchmark/`](benchmark/).

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
- 📧 **Email:** [wshuai@berkeley.edu](mailto:wshuai@berkeley.edu)  
- 🐛 **GitHub Issues:** open an issue in this repository

We’ll respond as soon as possible.

