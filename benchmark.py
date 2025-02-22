from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET



def read_ref(ref):
    REF = {}
    seq_dict = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos
    motif_sites = {}

    for r, contig in REF.items():
        print (r)
        for site in nt_search(str(contig), motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            motif_sites[r + ":" + str(site+exact_pos) + "+"] = motif_new

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            motif_sites[r + ":" + str(site+rev_exact_pos) + "-"] = motif_new

    return motif_sites

motif_new = "CTGCAG"
exact_pos = 5
ref = "/home/shuaiw/methylation/data/published_data/fanggang/bam/GCA_022869985.1_ASM2286998v1_genomic.fa"
REF = read_ref(ref)
motif_sites = get_motif_sites(REF, motif_new, exact_pos)
for site in motif_sites:
    print(site, motif_sites[site])

print (len(motif_sites))