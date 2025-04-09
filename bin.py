import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from typing import Dict, List
import os


def load_cluster_assignments(cluster_csv: str) -> Dict[int, List[str]]:
    """Load contig-cluster assignments from a CSV file."""
    cluster_df = pd.read_csv(cluster_csv)
    cluster_map = defaultdict(list)
    for _, row in cluster_df.iterrows():
        contig = str(row["contigs"])
        cluster = int(row["cluster"])
        cluster_map[cluster].append(contig)
    return cluster_map


def load_fasta_sequences(fasta_path: str) -> Dict[str, SeqIO.SeqRecord]:
    """Load all contig sequences into a dictionary."""
    return SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))


def write_bins_to_fasta(cluster_map: Dict[int, List[str]],
                        sequences: Dict[str, SeqIO.SeqRecord],
                        output_prefix: str = "bin") -> None:
    """Write each cluster/bin as a separate FASTA file."""
    for cluster_id, contigs in cluster_map.items():
        output_file = f"{output_prefix}_{cluster_id}.fasta"
        records = [sequences[ctg] for ctg in contigs if ctg in sequences]
        SeqIO.write(records, output_file, "fasta")
        print(f"✔️  Wrote {len(records)} sequences to {output_file}")


def bin_contigs_to_fastas(cluster_csv: str,
                          fasta_path: str,
                          output_prefix: str = "bin") -> None:
    """Main function to process input and create bin FASTA files."""
    ## check if the input files exist
    if not os.path.exists(cluster_csv):
        print (f"Cluster CSV file not found: {cluster_csv}")
        return
    cluster_map = load_cluster_assignments(cluster_csv)
    sequences = load_fasta_sequences(fasta_path)
    write_bins_to_fasta(cluster_map, sequences, output_prefix)
    print("✅ All bins written to FASTA files.")

if __name__ == "__main__":
    # Example usage:
    # cluster_csv = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_NM3/motif_cluster.h.csv"  # Path to the cluster CSV file
    # fasta_path = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa"  # Path to the input FASTA file
    # output_prefix = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_NM3/bins/bin"  # Prefix for output files

    cluster_csv = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_new/motif_cluster.csv"  # Path to the cluster CSV file
    fasta_path = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_200ppm_r1_LR_scaffold.fa"  # Path to the input FASTA file
    output_prefix = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_new/bins/bin"  # Prefix for output files

    bin_contigs_to_fastas(cluster_csv, fasta_path, output_prefix)

