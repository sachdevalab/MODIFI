"""
given a fasta file, use biopython to reverse one segment of it and output a new fasta file
"""

from Bio import SeqIO
from Bio.Seq import Seq


def reorder(input_fasta, output_fasta, breakpont):
    """
    Reverse complement a specific region in a FASTA file
    """
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            sequence = str(record.seq)
            
            # Extract the region to reverse complement
            before_region = sequence[:breakpont]
            after_region = sequence[breakpont:]
            
            # Reverse complement the region
            # Construct the new sequence
            new_sequence = after_region + before_region 
            record.id = record.id + "_reordered"
            # Create new record for the full sequence
            record.seq = Seq(new_sequence)
            SeqIO.write(record, outfile, "fasta")
            

if __name__ == "__main__":
    raw_fasta = "/home/shuaiw/borg/paper/E_faecalis/detail_3/infant_26_3_C.fa"
    new_fasta = "/home/shuaiw/borg/paper/E_faecalis/detail_3/circulize_infant_26_3_C.fa"
    reorder(raw_fasta, new_fasta, breakpont = 2436747)

    # raw_fasta = "/home/shuaiw/borg/paper/E_faecalis/detail_3/infant_14_31_C.fa"
    # new_fasta = "/home/shuaiw/borg/paper/E_faecalis/detail_3/circulize_infant_14_31_C.fa"
    # reorder(raw_fasta, new_fasta, breakpont = 2323856+84872)
