import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

### given a fasta file, which has several contigs
### break each contig into five parts, and output it in a new fasta file
### output the new fasta file


def break_contigs(fasta_file, break_fasta_file):
    output = open(break_fasta_file, "w")
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = record.seq
        seq_len = len(seq)
        if seq_len > 1000000:
            part_len = int(seq_len/5)
            for i in range(5):
                new_seq = seq[i*part_len:(i+1)*part_len]
                ## create a new empty record, not only copy the record
                new_id = record.id + "_" + str(i*part_len) + "_" + str((i+1)*part_len)
                new_record = SeqRecord(new_seq, id=new_id, description=record.description)
                # new_record.seq = new_seq
                # new_record.id = record.id + "_" + str(part_len) + "_" + str((i+1)*part_len)
                # new_record.description = ''# record.id + "_" + str(i)
                SeqIO.write(new_record, output, "fasta")
        else:
            SeqIO.write(record, output, "fasta")
    output.close()
    ## index the new fasta file with samtools
    cmd = "samtools faidx " + break_fasta_file
    
    os.system(cmd)
    return

## total contig number : 36608
# fasta_file = "/home/shuaiw/methylation/data/borg/contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs.fa"
# break_fasta_file = "/home/shuaiw/methylation/data/borg/contigs/all_break.contigs.fa"

fasta_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa"
break_fasta_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_break.fa"
break_contigs(fasta_file, break_fasta_file)