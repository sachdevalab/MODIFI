## given a assembly, select two contigs, break them into half, and generate two mosaic contigs

from Bio import SeqIO
import os

# Define file path and contig names
ref = "/home/shuaiw/borg/contigs/seven_contigs.fasta"
one = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_317_C"
second = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_961_C"

# Read sequences from FASTA file
contigs = {}
with open(ref, "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):

        contigs[record.id] = str(record.seq)

# Extract the selected contigs
seq_one = contigs.pop(one, None)  # Remove from dict
seq_second = contigs.pop(second, None)  # Remove from dict

# Ensure both contigs are found
if seq_one is None or seq_second is None:
    print("One or both contigs not found.")
    exit()

# Split each contig into two halves
half_len_one = len(seq_one) // 2
half_len_second = len(seq_second) // 2

first_half_one, second_half_one = seq_one[:half_len_one], seq_one[half_len_one:]
first_half_second, second_half_second = seq_second[:half_len_second], seq_second[half_len_second:]

# Generate mosaic contigs
mosaic_1 = first_half_one + second_half_second
mosaic_2 = first_half_second + second_half_one

# Output the sequences in FASTA format
output_file = "/home/shuaiw/borg/contigs/mis_contigs.fasta"
with open(output_file, "w") as out_fasta:
    # Write all original (unchanged) contigs
    for contig_id, seq in contigs.items():
        if contig_id == one or contig_id == second:
            print(f"skip {contig_id}")
        else:
            out_fasta.write(f">{contig_id}\n{seq}\n")
    
    # Write mosaic contigs
    out_fasta.write(f">{one}_mosaic\n{mosaic_1}\n")
    out_fasta.write(f">{second}_mosaic\n{mosaic_2}\n")

print(f"Modified contigs (including mosaics) saved to {output_file}")
## index the new contigs using samtools faidx
## samtools faidx mis_contigs.fasta
os.system(f"samtools faidx {output_file}")



