from Bio import SeqIO
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re




def read_bin_file(bin_file):
    """
    Read the bin file and return a dictionary of contig to bin.
    """
    bin_dict = {}
    with open(bin_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            contig = parts[0]
            bin_id = parts[1].strip()
            ## remove space 
            bin_id = bin_id.replace(" ", "")
            # if re.search("Ignavibacteriae_51_8_C", bin_id):
            #     print (f"Found Ignavibacteriae bin: {bin_id} for contig: {contig}")
            bin_dict[contig] = bin_id
    ## rewrite bin_file
    # with open(bin_file, 'w') as f:
    #     for contig, bin_id in bin_dict.items():
    #         f.write(f"{contig}\t{bin_id}\n")
    return bin_dict

def write_bin_fasta(assembly, bin_dict, bin_dir):
    ## use biopython to read the fasta file
    
    for record in SeqIO.parse(assembly, "fasta"):
        contig = record.id
        if contig in bin_dict:
            bin_id = bin_dict[contig]
            out_file = os.path.join(bin_dir, f"{bin_id}.fa")
            with open(out_file, 'a') as out_f:
                SeqIO.write(record, out_f, "fasta")


if __name__ == "__main__":
    assembly = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
    bin_file = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.bin.tab"
    bin_dir = "/home/shuaiw/methylation/data/borg/contigs/bins/"
    bin_dict = read_bin_file(bin_file)
    # write_bin_fasta(assembly, bin_dict, bin_dir)
    # print(f"Bin fasta files written to {bin_dir}")
