import os


def load_contigs():
    host_file = '/home/shuaiw/Methy/borg_test/borg.csv'
    ## read the contig name in a dict
    contig_dict = {}
    with open(host_file) as f:
        for line in f:
            contig = line.strip()
            contig_dict[contig] = 'borg'
    borg_file = '/home/shuaiw/Methy/borg_test/host.csv'
    with open(borg_file) as f:
        for line in f:
            contig = line.strip()
            contig_dict[contig] = 'host'
    control_dict = {}
    control_file = '/home/shuaiw/Methy/borg_test/control.csv'
    with open(control_file) as f:
        for line in f:
            contig = line.strip()
            if contig not in contig_dict:
                control_dict[contig] = 'control'

    contig_dict.update(control_dict)

    return contig_dict

## use biopython to extract the contigs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_contigs(whole_ref, contig_dict, extracted_ref):
    with open(extracted_ref, "w") as output_handle:
        for seq_record in SeqIO.parse(whole_ref, "fasta"):
            if seq_record.id in contig_dict:
                SeqIO.write(seq_record, output_handle, "fasta")
    ## index the fasta file
    os.system(f"samtools faidx {extracted_ref}")


contig_dict = load_contigs()
whole_ref = "/home/shuaiw/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
extracted_ref = "/home/shuaiw/borg/all_borgs2/all_borgs2.fa"
extract_contigs(whole_ref, contig_dict, extracted_ref)




# for contig in contig_dict:
#     motif_file = f"/home/shuaiw/methylation/data//borg/split_bam_dir2/{contig}.motif.csv"
#     motif_file2 = f"/home/shuaiw/methylation/data//borg/split_bam_dir/{contig}.motif.csv"
#     bam_file = f"/home/shuaiw/methylation/data//borg/split_bam_dir2/{contig}.bam"
#     fasta_file = f"/home/shuaiw/methylation/data//borg/split_bam_dir2/{contig}.fasta"
#     ## check if motif_file exists
#     # if os.path.exists(motif_file):
#     #     print (contig, motif_file)
#     #     # continue
#     # elif os.path.exists(motif_file2):
#     #     print (contig, motif_file2)
#     #     # continues
#     # else:
#     #     print (contig, "not exists")
#     #     # continue
#     ## check if bam exits
#     if os.path.exists(bam_file):
#         print (contig, bam_file)
#     if os.path.exists(fasta_file):
#         print (contig, fasta_file)
#     os.system(f'cp {bam_file}.pbi /home/shuaiw/methylation/data/borg/split_bam_dir3/')
#     # os.system(f'cp {fasta_file}.fai /home/shuaiw/methylation/data/borg/split_bam_dir3/')
    
