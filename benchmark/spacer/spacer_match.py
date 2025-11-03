# -*- coding: utf-8 -*-
import subprocess
import os
from pathlib import Path
import sys
import pandas as pd


def run_blastn_spacer_search(mge_fasta, outdir, spacer_fasta, raw_hit, min_id=0.95, threads=24):
    """
    Run BLASTN to find MGE targets of CRISPR spacers.
    
    Args:
        spacer_fasta (str): Path to spacer FASTA file (nucleotide).
        mge_fasta (str or Path): Path to MGE FASTA file (nucleotide).
        outdir (str or Path): Output directory to store results.
        min_id (float): Minimum percent identity (0-1). Default: 0.95.
        threads (int): Number of threads. Default: 24.
    
    Output:
        - Tabular alignment results in BLAST6 format.
    """
    mge_fasta = Path(mge_fasta).resolve()
    spacer_fasta = Path(spacer_fasta).resolve()
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    mgename = mge_fasta.stem
    hits_file = raw_hit
    blast_db_prefix = outdir / f"{mgename}_blastdb"

    # Step 1: Create BLAST database
    subprocess.run(
        f"makeblastdb  -in {mge_fasta} -dbtype nucl -out {blast_db_prefix}",
        shell=True,
        check=True
    )

    # Step 2: Run BLASTN search
    # subprocess.run(
    #     f"blastn -task blastn-short -query {spacer_fasta} -db {blast_db_prefix} -out {hits_file} "
    #     f"-evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' "
    #     f"-num_threads {threads} -perc_identity {min_id * 100:.1f}",
    #     shell=True,
    #     check=True
    # )

    ## rohan's version
    subprocess.run(
        f"blastn -task blastn-short -query {spacer_fasta} -db {blast_db_prefix} -out {hits_file} "
        f"-outfmt '6 std qlen nident' -max_target_seqs 5000 "
        f"-num_threads {threads}",
        shell=True,
        check=True
    )

    print(f"[✔] BLASTN search completed. Results written to: {hits_file}")

def predict_spacer(reference_fasta, outdir, prefix):
    """
    Predict CRISPR spacers in the reference genome using MinCED.
    
    Args:
        reference_fasta (str or Path): Path to reference genome FASTA file.
        spacer_fasta (str or Path): Path to output spacer FASTA file.
        outdir (str or Path): Output directory to store results.
    
    Output:
        - Spacer sequences in FASTA format.
    """
    reference_fasta = Path(reference_fasta).resolve()
    ## create output directory if not exist
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    # minced -spacers 96plex.hifiasm.p_ctg.rename.fa 96plex.hifiasm.p_ctg.rename.crisprs
    # Run MinCED to predict CRISPR spacers
    subprocess.run(
        f"minced -spacers {reference_fasta} {outdir}/{prefix}.crisprs",
        shell=True,
        check=True
    )
    ## print the running code
    # print(f"minced -spacers {reference_fasta} {outdir}/{prefix}.crisprs")

    print(f"[✔] Spacer prediction completed.")


def get_mge_fa(reference_fasta, mge_file, mge_fatsa):
    ## use biopython to extract the sequences from reference_fasta based on mge_file
    from Bio import SeqIO
    mge_ids = set()
    with open(mge_file, "r") as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            mge_id = fields[0]
            mge_ids.add(mge_id)
    ## read reference_fasta and write the sequences to mge_fatsa
    with open(reference_fasta, "r") as infile, open(mge_fatsa, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in mge_ids:
                SeqIO.write(record, outfile, "fasta")
    print(f"[✔] Extracted {len(mge_ids)} MGE sequences to {mge_fatsa}")
    return mge_ids

def filter_hit(raw_hit, spacer_linkage, mge_ids):
    """
    Filter BLAST hits to remove those where the spacer's contig is in the MGE set.
    Args:
        raw_hit (str): Path to raw BLAST hits file (tsv).
        spacer_linkage (str): Output file for filtered hits.
        mge_ids (set): Set of MGE contig IDs to exclude.
    Returns:
        pd.DataFrame: Filtered hits DataFrame.
    """
    import pandas as pd
    ## check if raw_hit is empty
    if not os.path.exists(raw_hit) or os.path.getsize(raw_hit) == 0:
        print(f"[⚠️] Raw hits file is empty or not found: {raw_hit}")
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    df = pd.read_csv(raw_hit, sep="\t", header=None)
    df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    
    # filter by length > 20 and pident > 95
    if "qseqid" not in df.columns:
        raise ValueError("Column 'qseqid' not found in BLAST hits file.")
    ## make sure the alignment length of qseqid is same as the length of qseqid
    new_df = pd.DataFrame()
    for index, row in df.iterrows():
        if row['length'] != row['qend'] - row['qstart'] + 1:
            new_df = new_df.append(row)
    df = new_df
    # Extract contig name from qseqid
    df['qseqid_ctg'] = df['qseqid'].str.split('_CRISPR_').str[0]
    # Remove hits where spacer contig is in MGE set
    df = df[~df["qseqid_ctg"].isin(mge_ids)]
    # Optionally filter by length and pident
    df = df[(df['length'] > 20) & (df['pident'] > 95)]
    # Save to spacer_linkage
    df.to_csv(spacer_linkage, sep="\t", index=False)
    print(f"[✔] Filtered hits saved to {spacer_linkage}, total {len(df)} hits.")
    return df


def filter_hit2(raw_hit, spacer_linkage, mge_ids):
    """
    Filter BLAST hits to remove those where the spacer's contig is in the MGE set.
    Args:
        raw_hit (str): Path to raw BLAST hits file (tsv).
        spacer_linkage (str): Output file for filtered hits.
        mge_ids (set): Set of MGE contig IDs to exclude.
    Returns:
        pd.DataFrame: Filtered hits DataFrame.
    """
    import pandas as pd
    ## check if raw_hit is empty
    if not os.path.exists(raw_hit) or os.path.getsize(raw_hit) == 0:
        print(f"[⚠️] Raw hits file is empty or not found: {raw_hit}")
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    df = pd.read_csv(raw_hit, sep="\t")
    # df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    df['qseqid_ctg'] = df['query_contig_id']
    # Remove hits where spacer contig is in MGE set
    df = df[~df["qseqid_ctg"].isin(mge_ids)]


    return df

def single_run():
    # prefix = sys.argv[1]
    # outdir = sys.argv[2]
    prefix = "ERR10042290"
    outdir = f"/groups/banfield/projects/multienv/methylation_temp/batch2_results/{prefix}"

    reference_fasta = f"{outdir}/{prefix}.hifiasm.p_ctg.rename.fa"
    spacer_outdir = f"{outdir}/spacer/"
    mge_file = f"{outdir}/all_mge.tsv"
    mge_fatsa = f"{outdir}/{prefix}_mge.fa"
    raw_hit = f"{spacer_outdir}/{Path(mge_fatsa).stem}_spacer_hits.tsv"
    filter_hit_file = f"{spacer_outdir}/{Path(mge_fatsa).stem}_spacer_hits.filter.tsv"
    spacer_linkage = f"{spacer_outdir}/{prefix}_spacer_linkage.tsv"
    # spacer_fasta = "/home/shuaiw/borg/paper/run2/all_spacer.fa"

    spacer_fasta = f"{spacer_outdir}/{prefix}_spacers.fa"
    predict_spacer(reference_fasta, spacer_outdir, prefix)
    
    mge_ids = get_mge_fa(reference_fasta, mge_file, mge_fatsa)
    run_blastn_spacer_search(mge_fasta=mge_fatsa, outdir=spacer_outdir, spacer_fasta=spacer_fasta, raw_hit=raw_hit, min_id=0.95, threads=24)
    cmd = f"""
    python3 /home/rohan/scripts/crispr_filter_blast.py -m 5 -s {spacer_fasta} -i {raw_hit} -o {filter_hit_file}
    """
    os.system(cmd)
    print (filter_hit_file)
    os.system(f"cat {filter_hit_file}")
    # filter_hit2(filter_hit_file, spacer_linkage, mge_ids)

def population_run():
    prefix = "population"
    resultdir = f"/groups/banfield/projects/multienv/methylation_temp/batch2_results/"
    outdir = f"/groups/banfield/projects/multienv/methylation_temp/batch2_spacer/"

    

    reference_fasta = f"{outdir}/{prefix}.all.fa"
    spacer_outdir = f"{outdir}/spacer/"
    mge_file = f"{outdir}/all_mge.tsv"
    mge_fatsa = f"{outdir}/{prefix}_mge.fa"
    raw_hit = f"{spacer_outdir}/{Path(mge_fatsa).stem}_spacer_hits.tsv"
    filter_hit_file = f"{spacer_outdir}/{Path(mge_fatsa).stem}_spacer_hits.filter.tsv"
    spacer_linkage = f"{spacer_outdir}/{prefix}_spacer_linkage.tsv"
    os.system(f"cat {resultdir}/*/*.hifiasm.p_ctg.rename.fa > {reference_fasta}")
    os.system(f"cat {resultdir}/*/all_mge.tsv > {mge_file}")

    spacer_fasta = f"{spacer_outdir}/{prefix}_spacers.fa"
    predict_spacer(reference_fasta, spacer_outdir, prefix)
    
    mge_ids = get_mge_fa(reference_fasta, mge_file, mge_fatsa)
    run_blastn_spacer_search(mge_fasta=mge_fatsa, outdir=spacer_outdir, spacer_fasta=spacer_fasta, raw_hit=raw_hit, min_id=0.95, threads=24)
    cmd = f"""
    python3 /home/rohan/scripts/crispr_filter_blast.py -m 5 -s {spacer_fasta} -i {raw_hit} -o {filter_hit_file}
    """
    os.system(cmd)
    # print (filter_hit_file)
    # os.system(f"cat {filter_hit_file}")
    df = filter_hit2(filter_hit_file, spacer_linkage, mge_ids)
    ## add taxonomy information to df
    df["query_taxa"] = df["query_contig_id"].apply(lambda x: get_taxa(x, isolation_taxa))
    df["target_taxa"] = df["target_id"].apply(lambda x: get_taxa(x, isolation_taxa))
    print (df)
    df.to_csv(spacer_linkage, sep="\t", index=False)
    os.system(f"cat {spacer_linkage}")
    print (spacer_linkage)

def get_taxa(contig, isolation_taxa):
    sra_id = contig.split("_")[0]
    if sra_id in isolation_taxa:
        return isolation_taxa[sra_id]
    else:
        return "unknown"

def read_gtdb(isolation_taxa, gatk):
    """
    Read the GTDB summary file and return a dictionary of contig to bin mapping.
    """
    gtdb_df = pd.read_csv(gatk, sep='\t')
    
    for index, row in gtdb_df.iterrows():
        anno = row['classification']
        sra_id = row['user_genome'].split("_")[0]
        isolation_taxa[sra_id] = anno
        return isolation_taxa
        
def load_all_taxa():
    isolation_taxa = {}
    resultdir = f"/groups/banfield/projects/multienv/methylation_temp/batch2_results/"
    for folder in os.listdir(resultdir):
        gatk = f"{resultdir}/{folder}/GTDB_2/gtdbtk.bac120.summary.tsv"
        print (gatk)
        ## skip if gatk file not exist
        if not os.path.exists(gatk):
            continue
        isolation_taxa = read_gtdb(isolation_taxa, gatk)
    # print (isolation_taxa)
    return isolation_taxa

if __name__ == "__main__":
    
    isolation_taxa = load_all_taxa()
    population_run()