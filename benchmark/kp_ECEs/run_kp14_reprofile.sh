#!/bin/bash
#SBATCH -p standard
#SBATCH -J kp14_reprofile
#SBATCH -o /home/shuaiw/borg/revision/kp_eces/gene_profile/reprofile_%j.out
#SBATCH --qos unlimit-submit-20-run
# Rebuild the 14-representative gene profile after restricting to the STRICT ECE set
# (kp14_pick_reps.py now filters members to ece_profile_final_sourcedata_ece.csv). Only the
# infant_2_845_L cluster representative changed (infant_4_361_C -> infant_2_845_L), but we regenerate
# the whole chain from kp14_representatives.tsv for consistency.
# Steps: gather rep nucl+ORFs -> eggNOG -> mob_typer/AMRFinder/abricate -> PLSDB ANI -> matrices.
set -uo pipefail
source /shared/software/miniconda3/latest/etc/profile.d/conda.sh 2>/dev/null || \
  source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
export PATH=/shared/software/bin:$PATH
T=${SLURM_CPUS_ON_NODE:-32}
GP=/home/shuaiw/borg/revision/kp_eces/gene_profile
SEQS=/home/shuaiw/borg/revision/kp_eces/seqs
RUN2=/groups/banfield/projects/multienv/methylation/data/borg/paper/run2
DBROOT=/home/shuaiw/borg/revision/kp_eces/dbs
SEQKIT="conda run -n tldr seqkit"
PRODIGAL=/shared/software/bin/prodigal
EMAPPER=/home/pengfanz/software/miniconda3/bin/emapper.py
ANICALC=/groups/sachdeva/projects/sag/SAGLink/workflow/anicalc.py
mkdir -p "$GP"
cd "$GP"

REPS=$(awk -F'\t' 'NR>1{print $2}' "$GP/kp14_representatives.tsv")
echo "[1] gather nucleotide for $(echo "$REPS"|wc -w) representatives"
: > kp14_reps.fna
for r in $REPS; do
  if $SEQKIT grep -p "$r" "$SEQS/all_metagenome_eces.fna" 2>/dev/null | grep -q .; then
    $SEQKIT grep -p "$r" "$SEQS/all_metagenome_eces.fna" >> kp14_reps.fna
  else
    samp=$(echo "$r" | sed -E 's/_[0-9]+_[CL]$//')
    $SEQKIT grep -p "$r" "$RUN2/$samp/$samp.hifiasm.p_ctg.rename.fa" >> kp14_reps.fna
  fi
done
echo "  reps in fna: $($SEQKIT stats -T kp14_reps.fna | awk 'NR==2{print $4}')"

echo "[2] predict ORFs (prodigal meta)"
$PRODIGAL -i kp14_reps.fna -a kp14_reps.faa -p meta -q -o /dev/null 2>/dev/null
# prodigal appends ' # start # end # ...' to headers; keep only the ORF id (<contig>_<n>)
sed -i 's/ .*//' kp14_reps.faa
echo "  ORFs: $(grep -c '^>' kp14_reps.faa)"

echo "[3] eggNOG-mapper (diamond)"
$EMAPPER -i kp14_reps.faa --itype proteins -m diamond \
  --data_dir /groups/diamond/databases/eggnog \
  --dmnd_db /groups/diamond/databases/eggnog/eggnog_proteins.dmnd \
  -o kp14 --output_dir "$GP" --cpu $T --override > eggnog.log 2>&1
echo "  annotated ORFs: $(grep -vc '^#' kp14.emapper.annotations 2>/dev/null)"

echo "[4] mob_typer / AMRFinder / abricate on representatives"
conda run -n kp_eces mob_typer --multi --infile kp14_reps.fna --out_file kp14_mobtyper.tsv \
  --num_threads $T --database_directory "$DBROOT/mob_suite" > mobtyper.log 2>&1
conda run -n kp_eces amrfinder -n kp14_reps.fna --plus --threads $T -o kp14_amrfinder.tsv > amrfinder.log 2>&1
/shared/software/bin/abricate --db vfdb kp14_reps.fna > kp14_abricate_vfdb.tsv 2> abricate.log
echo "  mob rows $(($(wc -l < kp14_mobtyper.tsv)-1)); amr rows $(($(wc -l < kp14_amrfinder.tsv)-1))"

echo "[5] PLSDB ANI for plasmid reps (viral reps keep existing viral ANI)"
conda run -n mod blastn -query kp14_reps.fna -db "$DBROOT/plsdb/plsdb" \
  -outfmt '6 std qlen slen' -max_target_seqs 25 -perc_identity 90 -num_threads 32 \
  > kp14_plsdb_blast.tsv 2> plsdb_blast.log
conda run -n methy3 python "$ANICALC" -i kp14_plsdb_blast.tsv -o kp14_plsdb_ani.tsv 2>> plsdb_blast.log
echo "  plsdb ani rows: $(($(wc -l < kp14_plsdb_ani.tsv)-1))"

echo "[6] matrices"
conda run -n methy3 python /home/shuaiw/MODIFI/benchmark/kp_ECEs/kp14_gene_matrix.py 2>&1 | tail -3
conda run -n methy3 python /home/shuaiw/MODIFI/benchmark/kp_ECEs/kp14_pfam_matrix.py 2>&1 | tail -2
echo "DONE reprofile"
