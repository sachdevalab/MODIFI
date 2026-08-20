#!/bin/bash
# Run Nanomotif MTase-linker on ONE isolate to get the reference-derived MTase table
# (mod_type_pred = ac[6mA/4mC] / m[5mC], motif_pred = predicted recognition sequence).
# Type-agnostic: we only consume the reference MTase predictions; the bin_motifs we feed
# is a minimal valid placeholder (mod_type='a') and its content does NOT drive the MTase typing.
#
# Usage: run_isolate_mtl.sh <ACC> <THREADS>
set -euo pipefail
ACC="$1"; THREADS="${2:-4}"
BATCH=/home/shuaiw/borg/paper/isolation/batch2_results
DEPS=/home/shuaiw/borg/revision/motif_class/nanomotif_mtl/ML_dependencies
OUTROOT=/home/shuaiw/borg/revision/motif_class/isolate_mtl
WORK="$OUTROOT/$ACC"

ASM="$BATCH/$ACC/$ACC.hifiasm.p_ctg.rename.fa"
MOTIFS="$BATCH/$ACC/${ACC}_methylation4/all.motifs.csv"
[ -s "$ASM" ] || { echo "[$ACC] no assembly"; exit 2; }
[ -s "$MOTIFS" ] || { echo "[$ACC] no motifs"; exit 3; }

# skip if already finished
if [ -s "$WORK/output/mtase_assignment_table.tsv" ]; then echo "[$ACC] already done, skip"; exit 0; fi

mkdir -p "$WORK"
# contig_bin.tsv: every contig -> bin=ACC
grep '^>' "$ASM" | sed 's/^>//; s/[[:space:]].*//' | awk -v b="$ACC" '{print $1"\t"b}' > "$WORK/contig_bin.tsv"

# bin_motifs.tsv: minimal valid placeholder from detected motifs (fraction>=0.5 & nDetected>=100).
python - "$MOTIFS" "$ACC" "$WORK/bin_motifs.tsv" <<'PY'
import sys, pandas as pd
from Bio.Seq import Seq
mf, acc, out = sys.argv[1], sys.argv[2], sys.argv[3]
df = pd.read_csv(mf)
df = df[(df.fraction>=0.5) & (df.nDetected>=100)]
rows=[]
seen=set()
for _,r in df.iterrows():
    m=str(r.motifString); p=int(r.centerPos)-1
    key=(m,p)
    if key in seen: continue
    seen.add(key)
    rc=str(Seq(m).reverse_complement())
    rows.append(dict(reference=acc, motif=m, mod_position=p, mod_type='a',
                     n_mod=100, n_nomod=1,
                     motif_type='non-palindrome', motif_complement=rc,
                     mod_position_complement=len(m)-1-p, n_mod_complement=100, n_nomod_complement=1))
cols=['reference','motif','mod_position','mod_type','n_mod','n_nomod','motif_type',
      'motif_complement','mod_position_complement','n_mod_complement','n_nomod_complement']
pd.DataFrame(rows, columns=cols).to_csv(out, sep='\t', index=False)
print(f"[{acc}] bin_motifs rows: {len(rows)}")
PY

nanomotif MTase-linker run \
  --assembly "$ASM" \
  --contig_bin "$WORK/contig_bin.tsv" \
  --bin_motifs "$WORK/bin_motifs.tsv" \
  -d "$DEPS" \
  -o "$WORK/output" \
  -t "$THREADS" > "$WORK/mtl.log" 2>&1

echo "[$ACC] done -> $WORK/output/mtase_assignment_table.tsv"
