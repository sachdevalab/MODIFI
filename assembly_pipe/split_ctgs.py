"""
Given a fasta file and a folder, split each contig to a fasta file in the folder
"""

from Bio import SeqIO
import os
import sys



def split_ctgs(fasta_file, out_folder):
    """
    Split each contig in the fasta file to a separate fasta file in the out_folder
    """

    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_id = record.id
        output_file = os.path.join(out_folder, f"{contig_id}.fasta")
        SeqIO.write(record, output_file, "fasta")
        # print(f"Written {contig_id} to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python split_ctgs.py <fasta_file> <out_folder>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    out_folder = sys.argv[2]

    split_ctgs(fasta_file, out_folder)
    print(f"Contigs from {fasta_file} have been split into {out_folder}")