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

def get_motif_sites(REF, motif_new, exact_pos, modified_loci, data):
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
                data.append([r, site, site+100, "red", motif_new, modified_loci[tag], "+"])
            else:
                data.append([r, site, site+100, "blue", motif_new, 0, "+"])
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
            motif_loci_num, motif_modify_num, ratio, proportion_all_modified], record_modified_sites, data



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
        if score < score_cutoff:
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


def get_df():
    REF = read_ref(my_ref)
    # print (REF)
    modified_loci = get_modified_ratio(gff)
    # motifs = pd.read_csv(all_motifs)
    ipd_ratio_dict = read_ipd_ratio(ipd_ratio_file)
    # print (len(ipd_ratio_dict))
    # data = []
    # for index, motif in motifs.iterrows():
        # motif_new = motif["motifString"]
        # exact_pos = motif["centerPos"]
    # motif_new = "AGCANNNNNNCCT"
    # exact_pos = 4
    motif_new = "GATC"
    exact_pos = 2
    all_record = {}
    data = []
    motif_profile, record_modified_sites, data = get_motif_sites(REF, motif_new, exact_pos, modified_loci, data)


    motif_new = "CTTCAG"
    exact_pos = 5
    motif_profile, record_modified_sites, data = get_motif_sites(REF, motif_new, exact_pos, modified_loci, data)

    motif_new = "AGCANNNNNNCCT"
    exact_pos = 4
    motif_profile, record_modified_sites, data = get_motif_sites(REF, motif_new, exact_pos, modified_loci, data)

    motif_new = "CAAYNNNNNCTGC"
    exact_pos = 3
    motif_profile, record_modified_sites, data = get_motif_sites(REF, motif_new, exact_pos, modified_loci, data)

    df = pd.DataFrame(data, columns=["chr", "start", "end", "color", "motif", "score", "strand"])
    df.to_csv("motif_df.csv", index=False)

def split_gff(gff2, contig, my_ref, out_dir):
    # bin_size = 1000 ##650  ## 50000
    bin_size = 55
    max_length = 50000000000
    name_dict = {}
    new_ref = out_dir + contig + ".fa"
    ## use seq IO to read the reference file
    for record in SeqIO.parse(my_ref, "fasta"):
        if len(record.seq) > max_length:
            record.seq = record.seq[:max_length]
        genome_length = len(record.seq)
        ## write the record to the new reference file
        SeqIO.write(record, new_ref, "fasta")
        break
    bin_size = genome_length // 1000 + 1
    # motif_list = ["CACAG", "GGAACG", "TACACG"]
    for motif in motif_list:
        # motif_gff = open(out_dir + contig + "_" + motif + ".gff", "w")
        motif_csv = out_dir + contig + "_" + motif + ".csv"
        data = []
        for line in open(gff2):
            if line[0] == "#":
                # motif_gff.write(line)
                continue
            line = line.strip().split("\t")
            if motif in line[8]:
                # motif_gff.write(line[0] + "\t" + line[1] + "\t" + motif + "\t" + line[3] + "\t" + line[4] + "\t" + line[5] + "\t" + line[6] + "\t" + line[7] + "\t" + line[8] + "\n")
                contig,type,start,stop,strand,score = line[0],motif,line[3],line[4],line[6],int(line[5])
                if int(start) >= max_length:
                    continue
                # score = score/60
                if score_cutoff >= score_cutoff:
                    score = 1
                else:
                    score = 0
                ## convert start and stop to bin
                bin_index = int(start) // bin_size
                start = round((bin_index + 0.25) * bin_size)
                stop = round((bin_index + 0.75) * bin_size)

                name = contig + ":" + str(start) + "-" + str(stop) + strand + "_" + motif
                if name not in name_dict:
                    data.append([name,contig,type,start,stop,strand, score, "yes"])
                    name_dict[name] = 1
        # motif_gff.close()
        df = pd.DataFrame(data, columns=["name", "contig", "type", "start", "stop", "strand", "score", "legend"])
        df.to_csv(motif_csv, index=False)


def all_split_gff(out_dir):
    bin_size = 5000 ##650  ## 50000
    max_length = 240707   ## E coli
    # max_length = 399533  ## K pneumoniae
    # max_length = 500000000  ## merged
    name_dict = {}
    # strain = "K_pneumoniae"
    # strain = "merged"
    strain = "Ecoli_H10407"
    new_ref = out_dir  + f"{strain}.fa"
    
    ## remove new_ref if it exists
    import os
    if os.path.exists(new_ref):
        os.remove(new_ref)
    
    # Open new_ref once for writing
    with open(new_ref, "w") as f_out:
        for contig in contig_list:
            my_ref = f"/home/shuaiw/borg/bench/zymo_new_ref/contigs/{contig}.fa"
            
            # Read and optionally truncate record, then write
            for record in SeqIO.parse(my_ref, "fasta"):
                if len(record.seq) > max_length:
                    record.seq = record.seq[:max_length]
                SeqIO.write(record, f_out, "fasta")
                break  # if only writing the first record from each file
    
    # motif_list = ["CACAG", "GGAACG", "TACACG"]
   
    for motif in motif_list:
        data = []
        new_csv = out_dir + f"{strain}_{motif}.csv"
        for contig in contig_list:
            gff2 = f"/home/shuaiw/borg/bench/zymo_new_ref/gffs/{contig}.reprocess.gff"
        
            
            for line in open(gff2):
                if line[0] == "#":
                    # motif_gff.write(line)
                    continue
                line = line.strip().split("\t")
                if motif in line[8]:
                    
                    contig,type,start,stop,strand,score = line[0],motif,line[3],line[4],line[6],int(line[5])
                    if int(start) >= max_length:
                        continue
                    if score < 30:
                        continue
                    ## convert start and stop to bin
                    # bin_index = int(start) // bin_size
                    # start = round((bin_index + 0.25) * bin_size)
                    # stop = round((bin_index + 0.75) * bin_size)

                    name = contig + ":" + str(start) + "-" + str(stop) + strand + "_" + motif
                    if name not in name_dict:
                        data.append([name,contig,type,start,stop,strand])
                        name_dict[name] = 1

        df = pd.DataFrame(data, columns=["name", "contig", "type", "start", "stop", "strand"])
        df.to_csv(new_csv, index=False)

if __name__ == "__main__":
    score_cutoff = 30
    contig = "infant_2_3_C"
    work_dir = f"/home/shuaiw/borg/paper/run2/infant_2/infant_2_methylation3/"
    motif_list = ["ATGCAT", "CAANNNNNNRTGA", "CAYNNNNNNTAYG"]
    my_ref = f"{work_dir}/contigs/{contig}.fa"
    gff2 = f"{work_dir}/gffs/{contig}.reprocess.gff"
    ipd_ratio_file = f"{work_dir}/ipd_ratio/{contig}.ipd3.csv"


    out_dir = "/home/shuaiw/borg/paper/circos/inversion/"
    split_gff(gff2, contig, my_ref, out_dir)


    # contig = "E_coli_H10407_6"
    # # contig = "B_cepacia_UCB-717_4"
    # motif_list = ["GATC", "CTTCAG"]
    # # motif_list = ["GATC", "CTTCAG", "AGCANNNNNNCCT", "CAAYNNNNNCTGC"]
    # my_ref = f"/home/shuaiw/borg/bench/zymo_new_ref/contigs/{contig}.fa"
    # gff2 = f"/home/shuaiw/borg/bench/zymo_new_ref/gffs/{contig}.reprocess.gff"
    # ipd_ratio_file = f"/home/shuaiw/borg/bench/zymo_new_ref_NM3/ipd_ratio/{contig}.ipd3.csv"


    # motif_list = ["CAGAC", "CCGG", "TGCCCA", "TCTANNNNNNNRTNG","GAANNNNNNTGGC"]
    # # contig = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_553_L"
    # contig = "SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_10354_C"
    # my_ref = f"/home/shuaiw/borg/bench/soil/run1/contigs/{contig}.fa"
    # gff2 = f"/home/shuaiw/borg/bench/soil/run1/gffs/{contig}.reprocess.gff"
    # ipd_ratio_file = f"/home/shuaiw/borg/bench/soil/run1/ipd_ratio/{contig}.ipd3.csv"


    # out_dir = "/home/shuaiw/borg/paper/circos/inversion"
    # motif_list = ["GATC", "CTTCAG", "AGCANNNNNNCCT", "CAAYNNNNNCTGC"]
    # contig_list = ["E_coli_H10407_1", "E_coli_H10407_2", "E_coli_H10407_3","E_coli_H10407_4","E_coli_H10407_5","E_coli_H10407_6"]
    # motif_list = ["GATC", "CGCATC"]  ## K_pneumoniae_BAA-2146_1
    # contig_list = ["K_pneumoniae_BAA-2146_1", "K_pneumoniae_BAA-2146_2","K_pneumoniae_BAA-2146_3", "K_pneumoniae_BAA-2146_4", "K_pneumoniae_BAA-2146_5"]

    # motif_list = ["TTGANNNNNNCCT", "CGTCGVNY", "ACAYNNNNNNNTGNG"]  
    # contig_list = ["B_cereus_971_1","B_multivorans_249_1"]

    # split_gff(gff2, contig, my_ref, out_dir)
    # all_split_gff(out_dir)



    




