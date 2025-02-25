from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
# import matplotlib.pyplot as plt
import sys


def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        # return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos, modified_loci):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos + 1
    motif_sites = {}
    motif_loci_num = 0
    motif_modify_num = 0
    for_loci_num = 0
    rev_loci_num = 0
    for_modified_num = 0
    rev_modified_num = 0

    for r, contig in REF.items():
        for site in nt_search(str(contig), motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            tag = r + ":" + str(site+exact_pos) + "+"
            motif_loci_num += 1
            for_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                for_modified_num += 1

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            tag = r + ":" + str(site+rev_exact_pos) + "-"
            motif_loci_num += 1
            rev_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                rev_modified_num += 1
    # print ("for_loci_num", for_loci_num, for_modified_num, "forward modified ratio", for_modified_num/for_loci_num)
    # print ("rev_loci_num", rev_loci_num, rev_modified_num, "reverse modified ratio", rev_modified_num/rev_loci_num)

    # print ("motif_loci_num", motif_loci_num)
    # print ("motif_modify_num", motif_modify_num)
    if motif_loci_num == 0:
        ratio = 0
    else:
        ratio = motif_modify_num/motif_loci_num
    if for_loci_num == 0:
        for_ratio = 0
    else:
        for_ratio = for_modified_num/for_loci_num
    if rev_loci_num == 0:
        rev_ratio = 0
    else:
        rev_ratio = rev_modified_num/rev_loci_num
    if len(modified_loci) == 0:
        proportion_all_modified = 0
    else:
        proportion_all_modified = motif_modify_num/len(modified_loci)

    if ratio > 0.4 and proportion_all_modified > 0.1:
        print (motif_new, "modified_ratio", ratio, motif_modify_num, motif_loci_num, proportion_all_modified)

    return [for_loci_num, for_modified_num,for_ratio,\
            rev_loci_num, rev_modified_num, rev_ratio,\
            motif_loci_num, motif_modify_num, ratio, proportion_all_modified]



def get_modified_ratio(gff):
    ## read the gff file
    f = open(gff, "r")
    modified_loci = {}
    for line in f:
        if line[0] == "#":
            continue
        line = line.strip().split("\t")

        ref = line[0]
        pos = int(line[3]) 
        strand = line[6]
        score = int(line[5])
        if score <= score_cutoff:
            continue
        modified_loci[ref + ":" + str(pos) + strand] = score
    print ("no. of modified loci", len(modified_loci))
    return modified_loci

         
if __name__ == "__main__":
    # motif_new = "CTGCAG"
    # exact_pos = 5
    score_cutoff = 30
    # my_ref = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
    # gff = "/home/shuaiw/borg/bench/C227_native/gffs/CP011331.1.gff"

    # my_ref = "/home/shuaiw/borg/all_test//contigs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_19121_L.fa"
    # gff = "/home/shuaiw/borg/all_test//gffs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_19121_L.gff"
    # all_motifs = "/home/shuaiw/borg/all_test/test_motifs.csv"
    # profile = "/home/shuaiw/borg/all_test/profiles/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_19121_L.csv"

    my_ref = sys.argv[1]
    gff = sys.argv[2]
    all_motifs = sys.argv[3]
    profile = sys.argv[4]

    REF = read_ref(my_ref)
    # print (REF)
    modified_loci = get_modified_ratio(gff)
    motifs = pd.read_csv(all_motifs)
    data = []
    for index, motif in motifs.iterrows():
        motif_new = motif["motifString"]
        exact_pos = motif["centerPos"]
        motif_profile = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
        data.append([motif_new, exact_pos] + motif_profile)

    df = pd.DataFrame(data, columns = ["motifString", "centerPos", "for_loci_num", "for_modified_num", "for_modified_ratio",\
                                        "rev_loci_num", "rev_modified_num", "rev_modified_ratio",\
                                        "motif_loci_num", "motif_modified_num", "motif_modified_ratio", "proportion"])
    df.to_csv(profile, index=False)


