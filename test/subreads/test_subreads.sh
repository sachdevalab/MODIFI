#!/bin/bash
## Resolve inputs, control_db and main.py relative to this script so it works both
## from a cloned repo and from the installed share/modifi/test/subreads/ location.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

## remove output directory (in the current working directory) if it exists
rm -rf output/

python "$SCRIPT_DIR/../../main.py" --aligned_bam "$SCRIPT_DIR/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_100_L.bam" \
-o output/ \
-r "$SCRIPT_DIR/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_100_L.fa" \
--read_type subreads \
--kmer_mean_db "$SCRIPT_DIR/../../control_db/control_db.up7.down3.mean.dat" \
--kmer_num_db "$SCRIPT_DIR/../../control_db/control_db.up7.down3.num.dat" \
--threads 5