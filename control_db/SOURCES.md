# Control database sources and provenance

MODIFI normalizes sample IPD against a per-k-mer (`up7.down3`, k=10, 4^10 = 1,048,576 entries)
control model. PacBio chemistries differ in baseline kinetics, so the control database should be
**chemistry-matched** to the sample. This file records what each shipped database is, where it
came from, and its k-mer coverage, so the provenance stays traceable (the binaries alone do not
encode chemistry).

Each database is a pair: `*.mean.dat` (per-k-mer mean IPD) and `*.num.dat` (per-k-mer observation
count), both 4,194,312 bytes. Identify a sample's chemistry from its BAM header:

```bash
samtools view -H input.bam | grep '@RG' | tr ';' '\n' | grep -E 'BINDINGKIT|SEQUENCINGKIT' | sort -u
```

| Chemistry | BINDINGKIT (example) | Files | Pass via |
|---|---|---|---|
| Sequel II S/P5-C2 | 101-894-200 / 102-194-100 | `control_db.up7.down3.{mean,num}.dat` (default) | `--kmer_mean_db` / `--kmer_num_db` |
| Revio R/P1-C1 | 102-739-100 | `control_db.Revio.up7.down3.{mean,num}.dat` | `--kmer_mean_db` / `--kmer_num_db` |
| RS II (subreads) | 100236500 | `control_db.RSII.up7.down3.{mean,num}.dat` | `--kmer_mean_db` / `--kmer_num_db` |

---

## Sequel II S/P5-C2 (default, unlabeled files)

- **Files:** `control_db.up7.down3.mean.dat`, `control_db.up7.down3.num.dat`
- **Chemistry:** Sequel II, S/P5-C2 (HiFi). This is the default DB; Sequel II users can use it as-is.
- **Source:** 96plex Zymo mock community metagenome
  (`run2/96plex/96plex_methylation4/control/`).
- **Read type:** hifi
- **k-mer coverage:** 1,048,576 / 1,048,576 nonzero (100.0%); total observations 186,816,991;
  mean 178.2 obs per k-mer.
- **md5:** mean `37885cfe8706aad9b6cfab07e46741f0`, num `1d77bf1c013f3c0e7cc2c05f852c841f`
- **Note:** built from a low-complexity mock, so per-k-mer depth (~178) is far lower than a
  high-complexity metagenome. Adequate as a default, but a soil-derived Sequel II S/P5-C2 control
  (e.g. `run2/soil_2`, ~2400 obs/k-mer) is available if higher-depth statistics are needed.

## Revio R/P1-C1

- **Files:** `control_db.Revio.up7.down3.mean.dat`, `control_db.Revio.up7.down3.num.dat`
- **Chemistry:** Revio, R/P1-C1 (HiFi).
- **Source:** high-complexity soil metagenome
  (`run2/soil_s3_1/soil_s3_1_methylation4/control/`), copied verbatim.
- **Read type:** hifi
- **k-mer coverage:** 1,048,576 / 1,048,576 nonzero (100.0%); total observations 2,129,546,881;
  mean 2030.9 obs per k-mer.
- **md5:** mean `dc2f0290ad363251725f9cec6448b8c8`, num `dea60a80667cdb6e0f868c1a566ec24c`

## RS II (subreads)

- **Files:** `control_db.RSII.up7.down3.mean.dat`, `control_db.RSII.up7.down3.num.dat`
- **Chemistry:** PacBio RS II, subreads (modification-free WGA control).
- **Source:** 3-strain whole-genome-amplified (WGA) control, aligned subreads:
  - BAM: `published_data/fanggang/align/merge_WGA.align.bam`
  - Reference: `published_data/fanggang/ref/three_species.fa`, 3 contigs (~10.7 Mbp):
    - `CP011330.1` Helicobacter pylori J99 (1.70 Mbp)
    - `CP011331.1` Escherichia coli O104:H4 str. C227-11 (5.29 Mbp)
    - `NC_007963.1` Chromohalobacter salexigens DSM 3043 (3.70 Mbp)
- **Read type:** subreads
- **k-mer coverage:** 1,030,306 / 1,048,576 nonzero (98.3%); total observations 16,634,795;
  mean 16.1 obs per k-mer. Coverage/depth is lower than the metagenome DBs (inherent to a
  3-isolate control), but nearly all 10-mers are represented and the WGA input is
  modification-free.
- **md5:** mean `ce9eb5d26001222ea48da40109178945`, num `b5afdac3c4367eba79b9eaf0278935cf`
- **Build (2026-09-07), MODIFI_subreads env:**
  ```bash
  python main.py \
    --aligned_bam published_data/fanggang/align/merge_WGA.align.bam \
    -r published_data/fanggang/ref/three_species.fa \
    -o <workdir> --read_type subreads --up 7 --down 3 --threads 48
  # then copy <workdir>/control/control_db.up7.down3.{mean,num}.dat here as control_db.RSII.*
  ```
- **Consumed by:** `benchmark/motif_classification/jf8_wga.sh` (chemistry-matched control for the
  JF8 mock).

---

### Notes

- The k-mer window (`up7.down3`) must match between sample and control.
- **Rebuilding a control DB:** MODIFI's end-of-run cleanup deletes the run's own `<output>/control/`
  directory (`main.py`, cleanup step) when a full pipeline completes without `--no-clean`. When you
  build a reusable control DB, copy `control/control_db.up7.down3.{mean,num}.dat` out **before** the
  run finishes, or pass `--no-clean`, otherwise the DB is removed with the other intermediates.
