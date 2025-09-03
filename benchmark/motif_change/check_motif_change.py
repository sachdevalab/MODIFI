import profile
from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
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

    motif_ipd_ratio = []
    record_modified_sites = {}

    for r, contig in REF.items():
        for site in nt_search(str(contig), motif_new)[1:]:
            # for i in range(site, site + motif_len):
            #     motif_sites[r + ":" + str(i) + "+"] = motif_new
            tag = r + ":" + str(site+exact_pos) + "+"
            record_modified_sites[tag] = motif_new
            ##  also mark the nearby sites
            # for i in range(site-100, site+100):
            #     record_modified_sites[r + ":" + str(i) + "+"] = motif_new
            # print (tag)
            motif_loci_num += 1
            for_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                for_modified_num += 1
            if tag in ipd_ratio_dict:
                motif_ipd_ratio.append(ipd_ratio_dict[tag][1])

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            tag = r + ":" + str(site+rev_exact_pos) + "-"
            record_modified_sites[tag] = motif_new
            # for i in range(site-100, site+100):
            #     record_modified_sites[r + ":" + str(i) + "+"] = motif_new
            motif_loci_num += 1
            rev_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                rev_modified_num += 1
            if tag in ipd_ratio_dict:
                motif_ipd_ratio.append(ipd_ratio_dict[tag][1])
    # print ("for_loci_num", for_loci_num, for_modified_num, "forward modified ratio", for_modified_num/for_loci_num)
    # print ("rev_loci_num", rev_loci_num, rev_modified_num, "reverse modified ratio", rev_modified_num/rev_loci_num)

    ## plot the distribution of motif_ipd_ratio
    print (len(motif_ipd_ratio))
    # plt.hist(motif_ipd_ratio, bins=100)
    # plt.xlabel("IPD ratio")
    # plt.ylabel("Frequency")
    # plt.title("IPD ratio distribution of motif " + motif_new)
    # plt.savefig("../tmp/" + motif_new + ".png")
    # plt.close()

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
            motif_loci_num, motif_modify_num, ratio, proportion_all_modified], record_modified_sites

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

def read_ipd_ratio(ipd_ratio_file):
    ipd_ratio_dict = {}
    df = pd.read_csv(ipd_ratio_file, nrows = 1000)
    for index, row in df.iterrows():
        # print (row['refName'], row['tpl'], row['strand'], row['coverage'], row['ipd_ratio'])
        if row['strand'] == 1:
            strand_string = "-"
        else:
            strand_string = "+"
        tag = row['refName'] + ":" + str(row['tpl']+1) + strand_string
        ipd_ratio_dict[tag] = [row['ipd_ratio'], row['tMean']]
    print (len(ipd_ratio_dict))
    ## print some of the ipd_ratio_dict
    for i, (k, v) in enumerate(ipd_ratio_dict.items()):
        # print (k, v)
        if i > 10:
            break
    return ipd_ratio_dict

         
if __name__ == "__main__":
    # motif_new = "CTGCAG"
    # exact_pos = 5
    score_cutoff = 30
    # my_ref = "/home/shuaiw/methylation/data/published_data/fanggang/ref/C227.fa"
    # gff = "/home/shuaiw/borg/bench/C227_native/gffs/CP011331.1.gff"

    # my_ref = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/contigs/E_coli_H10407_1.fa"
    # gff = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/gffs/E_coli_H10407_1.gff"
    # # gff = "/home/shuaiw/borg/bench/test/E_coli_H10407_1.gff"
    # # all_motifs = "/home/shuaiw/methylation/data/borg/bench/zymo2/all.motifs.csv"
    # # profile = "/home/shuaiw/borg/test.csv"
    # ipd_ratio_file = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30/ipd_ratio/E_coli_H10407_1.ipd3.csv"

    prefix_list = [["cow_bioreactor_1", "cow_bioreactor_1_636_C"], \
                   ["cow_bioreactor_2", "cow_bioreactor_2_1062_L"], \
                   ["cow_bioreactor_2", "cow_bioreactor_2_601_C"], \
                   ["cow_bioreactor_4", "cow_bioreactor_4_1750_C"],\
                    ["cow_bioreactor_5", "cow_bioreactor_5_1162_C"]]
    # prefix = "cow_bioreactor_4"
    # contig = "cow_bioreactor_4_1750_C"
    data = []
    motif_list = [["GATC", 2], ["ACNCAG", 5], ["GAAATC", 4], ["ACTNNNNNNRGTC", 1], ["GGCATC", 4]]
    for prefix, contig in prefix_list:
        my_ref = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3/contigs/{contig}.fa"
        gff = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3/gffs/{contig}.gff"
        ipd_ratio_file = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation3/ipd_ratio/{contig}.ipd3.csv"
        REF = read_ref(my_ref)
        # print (REF)
        modified_loci = get_modified_ratio(gff)
        # motifs = pd.read_csv(all_motifs)
        ipd_ratio_dict = read_ipd_ratio(ipd_ratio_file)
        
        for motif_new, exact_pos in motif_list:
            all_record = {}
            motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
            print (motif_profile)
            data.append([contig,motif_new, motif_profile[-2]])
    df = pd.DataFrame(data, columns = ["contig", "motifString", "fraction"])
    ## plot heatmap with df
    df_pivot = df.pivot(index='contig', columns='motifString', values='fraction').fillna(0)
    import seaborn as sns
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(df_pivot, cmap="YlGnBu", annot=True, fmt=".2f", cbar_kws={'label': 'Fraction'})
    ax.set_title("Motif Fraction Heatmap")
    ax.set_xlabel("Motif String")
    ax.set_ylabel("Contig")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()
    plt.savefig("../../tmp/results2/motif_fraction_heatmap.pdf")

    # my_ref = sys.argv[1]
    # gff = sys.argv[2]
    # all_motifs = sys.argv[3]
    # profile = sys.argv[4]

    # REF = read_ref(my_ref)
    # # print (REF)
    # modified_loci = get_modified_ratio(gff)
    # # motifs = pd.read_csv(all_motifs)
    # ipd_ratio_dict = read_ipd_ratio(ipd_ratio_file)
    # # print (len(ipd_ratio_dict))
    # # data = []
    # # for index, motif in motifs.iterrows():
    #     # motif_new = motif["motifString"]
    #     # exact_pos = motif["centerPos"]
    # # motif_new = "AGCANNNNNNCCT"
    # # exact_pos = 4
    # # motif_new = "GGCATC"
    # # exact_pos = 4
    # motif_list = [["GATC", 2], ["ACNCAG", 5], ["GAAATC", 4], ["ACTNNNNNNRGTC", 1], ["GGCATC", 4]]
    # for motif_new, exact_pos in motif_list:
    #     all_record = {}
    #     motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
    #     print (motif_profile)
    # ## add record_modified_sites to all_record
    # for k, v in record_modified_sites.items():
    #     all_record[k] = v
    # motif_new = "CTTCAG"
    # exact_pos = 5
    # motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
    # ## add record_modified_sites to all_record
    # for k, v in record_modified_sites.items():
    #     all_record[k] = v
    # motif_new = "CTGAAG"
    # exact_pos = 5
    # motif_profile, record_modified_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)
    # ## add record_modified_sites to all_record
    # for k, v in record_modified_sites.items():
    #     all_record[k] = v
    # record_modified_sites  = all_record
    # print ("no. of final modified sites", len(record_modified_sites))
    
    # # print (motif_new, motif_profile)
    # new_gff = "/home/shuaiw/borg/bench/test/E_coli_H10407_1.new.gff"
    # h = open(new_gff, "w")
    # ## read the gff file
    # f = open(gff, "r")
    # modified_loci = {}
    # for line in f:
    #     if line[0] == "#":
    #         print (line.strip(), file = h)
    #         continue
    #     field = line.strip().split("\t")

    #     ref = field[0]
    #     pos = int(field[3]) 
    #     strand = field[6]
    #     score = int(field[5])
    #     tag = ref + ":" + str(pos) + strand
    #     if tag not in all_record and score >= 40:
    #         print (line.strip(), file = h)
    # h.close()


    ## filter gff file

    #     data.append([motif_new, exact_pos] + motif_profile)

    # df = pd.DataFrame(data, columns = ["motifString", "centerPos", "for_loci_num", "for_modified_num", "for_modified_ratio",\
    #                                     "rev_loci_num", "rev_modified_num", "rev_modified_ratio",\
    #                                     "motif_loci_num", "motif_modified_num", "motif_modified_ratio", "proportion"])
    # df.to_csv(profile, index=False)


