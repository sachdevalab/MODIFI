from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

def rotate_genome(input_fasta, output_fasta):
    # Read all contigs from the input FASTA
    records = list(SeqIO.parse(input_fasta, "fasta"))

    # Process the first contig
    first_record = records[0]
    seq = str(first_record.seq)
    length = len(seq)
    new_start = int(length / 2)

    if not (0 < new_start < length):
        raise ValueError(f"new_start must be between 1 and {length-1}")

    # Rotate the sequence of the first contig
    new_seq = seq[new_start:] + seq[:new_start]

    # Update the first record
    first_record.seq = new_seq
    first_record.id = first_record.id + f"_rotated_{new_start}"
    first_record.description = f"Rotated circular genome, new start at position {new_start}"

    # Save all contigs to the output FASTA
    with open(output_fasta, "w") as f:
        SeqIO.write(records, f, "fasta")

    print(f"Genome re-linearized and saved to {output_fasta}")

# Example usage:
# input_fasta = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/asems/E_coli_H10407_bc2007.fa"
# output_fasta = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/re_linearize/E_coli_H10407_bc2007_re_linearized.fa"
input_fasta = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/asems/E_coli_K12-MG1655_bc2008.fa"
output_fasta = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/re_linearize/E_coli_K12-MG1655_bc2008_re_linearized.fa"
rotate_genome(input_fasta, output_fasta)