"""
given a fasta file, use biopython to reverse one segment of it and output a new fasta file
"""

from Bio import SeqIO
from Bio.Seq import Seq

raw_fasta = "/home/shuaiw/borg/paper/E_faecalis/infant_14_31_C.fa"
new_fasta = "/home/shuaiw/borg/paper/E_faecalis/infant_14_31_C_reverse.fa"
inversion_fasta = "/home/shuaiw/borg/paper/E_faecalis/inversion.fa"

# Define the region to reverse complement
start_pos = 2323857 - 1  # Convert to 0-based indexing
end_pos = 2323857 + 3259 - 1  # 2327116 in 0-based indexing

def reverse_complement_region(input_fasta, output_fasta, start, end, inversion_fasta):
    """
    Reverse complement a specific region in a FASTA file
    """
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile, open(inversion_fasta, "w") as inv_file:
        for record in SeqIO.parse(infile, "fasta"):
            sequence = str(record.seq)
            
            # Extract the region to reverse complement
            before_region = sequence[:start]
            region_to_reverse = sequence[start:end+1]
            after_region = sequence[end+1:]
            
            # Reverse complement the region
            reversed_region = str(Seq(region_to_reverse).reverse_complement())
            
            # Construct the new sequence
            new_sequence = before_region + reversed_region + after_region
            
            # Create new record for the full sequence
            record.seq = Seq(new_sequence)
            SeqIO.write(record, outfile, "fasta")
            
            # Create a separate record for just the inverted sequence
            inversion_record = SeqIO.SeqRecord(
                Seq(reversed_region),
                id=f"{record.id}_inversion_{start+1}_{end+1}",
                description=f"Inverted region from {start+1} to {end+1} of {record.id}"
            )
            SeqIO.write(inversion_record, inv_file, "fasta")
            
            print(f"Reversed region {start+1}-{end+1} (1-based coordinates)")
            print(f"Original region: {region_to_reverse[:50]}...")
            print(f"Reversed region: {reversed_region[:50]}...")
            print(f"Inversion sequence saved to: {inversion_fasta}")

if __name__ == "__main__":
    reverse_complement_region(raw_fasta, new_fasta, start_pos, end_pos, inversion_fasta)
