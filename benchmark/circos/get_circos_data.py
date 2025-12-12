from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys
import os
import profile
from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys, os, re
import numpy as np
from collections import defaultdict
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_unique_motifs, My_sample, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa




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


def get_all_loci(gff2, contig_name, my_ref, out_dir,motif_list, depth="NA"):
    modified_loci = get_modified_ratio(gff2)
    new_ref = os.path.join(out_dir, contig_name + ".fa")
    REF = read_ref(my_ref)

    ## use seq IO to read the reference file
    for record in SeqIO.parse(my_ref, "fasta"):

        genome_length = len(record.seq)
        ## write the record to the new reference file
        SeqIO.write(record, new_ref, "fasta")
        break
    bin_size = genome_length // 1000 + 1


    all_data = []
    for motif_data in motif_list:
        motif_new, exact_pos = motif_data
        motif_csv = os.path.join(out_dir, contig_name + "_" + motif_new + ".csv")
        data = []
        
        rev_exact_pos = len(motif_new) - exact_pos + 1
        for r, contig in REF.items():
            for site in nt_search(str(contig), motif_new)[1:]:
                position = site + exact_pos
                tag = r + ":" + str(position) + "+"
                name = contig_name + ":" + str(position) + "-" + str(position) + "+" + "_" + motif_new
                if tag in modified_loci:
                    data.append([name,contig_name,motif_new,position,position,"+", 1, "yes"])
                else:
                    data.append([name,contig_name,motif_new,position,position,"+", 0.5, "no"])


            for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
                1:
            ]:
                position = site + rev_exact_pos
                tag = r + ":" + str(position) + "-"
                name = contig_name + ":" + str(position) + "-" + str(position) + "-" + "_" + motif_new
                if tag in modified_loci:
                    data.append([name,contig_name,motif_new,position,position,"-", 1, "yes"])
                else:
                    data.append([name,contig_name,motif_new,position,position,"-", 0.5, "no"])

        all_data.extend(data)   
        df = pd.DataFrame(data, columns=["name", "contig", "type", "start", "stop", "strand", "score", "legend"])
        df.to_csv(motif_csv, index=False)
    all_df = pd.DataFrame(all_data, columns=["name", "contig", "type", "start", "stop", "strand", "score", "legend"])
    plot_line(all_df, out_dir, new_ref, contig_name, depth=depth)

def calculate_gc_skew(sequence, window_size=1000):
    """Calculate GC skew for a sequence with sliding window"""
    gc_skew = []
    positions = []
    
    for i in range(0, len(sequence) - window_size + 1, window_size // 2):
        window = sequence[i:i + window_size]
        g_count = window.count('G') + window.count('g')
        c_count = window.count('C') + window.count('c')
        
        if g_count + c_count > 0:
            skew = (g_count - c_count) / (g_count + c_count)
        else:
            skew = 0
            
        gc_skew.append(skew)
        positions.append(i + window_size // 2)
    
    return positions, gc_skew

def plot_line(df, out_dir, new_ref, contig_name, depth="NA"):
    ## for each motif, for each strand, plot a barplot
    ## the y is the score, save the plot
    import matplotlib.pyplot as plt
    
    sns.set(style="whitegrid")
    out_dir = out_dir + "/figs"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    with PdfPages(os.path.join(out_dir, contig_name + "_motif_plots.pdf")) as pdf:
        # Get genome length for consistent x-axis scaling
        genome_length = 0
        for record in SeqIO.parse(new_ref, "fasta"):
            genome_length = len(record.seq)
            positions, gc_skew = calculate_gc_skew(str(record.seq))
            
            plt.figure(figsize=(25, 2))
            
            # Convert to numpy arrays for easier manipulation
            positions_array = np.array(positions)
            gc_skew_array = np.array(gc_skew)
            
            # Fill positive values (above zero) in green
            plt.fill_between(positions_array, 0, gc_skew_array, 
                           where=(gc_skew_array >= 0), 
                           color='green', alpha=0.7, label='Positive GC Skew')
            
            # Fill negative values (below zero) in purple/magenta
            plt.fill_between(positions_array, 0, gc_skew_array, 
                           where=(gc_skew_array < 0), 
                           color='purple', alpha=0.7, label='Negative GC Skew')
            
            # Add zero line
            plt.axhline(y=0, color='black', linestyle='-', alpha=0.8, linewidth=0.5)
            
            plt.title(f"GC Skew - {record.id} (Depth: {depth})")
            plt.ylabel("GC Skew")
            plt.grid(True, alpha=0.3)
            plt.legend()
            
            # Set consistent x-axis limits and ticks
            plt.xlim(0, genome_length)
            ax = plt.gca()
            x_ticks = np.linspace(0, genome_length, 21)  # 20 intervals = 21 ticks
            ax.set_xticks(x_ticks)
            ax.set_xticklabels([f"{int(x):,}" for x in x_ticks], rotation=45)
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
            break  # Only process first sequence

        for strand in df['strand'].unique():
            for motif in df['type'].unique():
                subset = df[(df['type'] == motif) & (df['strand'] == strand)]
                print (strand, motif, len(subset))
                if len(subset) > 0:  # Only plot if there's data
                    plt.figure(figsize=(25, 2))
                    
                    # Use scatter plot with different colors for better visualization
                    modified = subset[subset['score'] == 1]
                    unmodified = subset[subset['score'] == 0.5]
                    
                    if len(unmodified) > 0:
                        plt.scatter(unmodified['start'], unmodified['score'], 
                                  color='blue', alpha=0.6, s=10, label='Unmodified')
                    if len(modified) > 0:
                        plt.scatter(modified['start'], modified['score'], 
                                  color='red', alpha=0.8, s=15, label='Modified')
                    
                    plt.title(f"Motif: {motif}, Strand: {strand} (n={len(subset)})")
                    plt.ylabel("Modification Score")
                    plt.ylim(0.3, 1.2)  # Add some padding around the discrete values
                    
                    if len(modified) > 0 and len(unmodified) > 0:
                        plt.legend()
                    
                    # Set same x-axis scale as GC skew plot
                    plt.xlim(0, genome_length)
                    # ax = plt.gca()
                    # x_ticks = np.linspace(0, genome_length, 21)  # 20 intervals = 21 ticks
                    # ax.set_xticks(x_ticks)
                    # ax.set_xticklabels([f"{int(x):,}" for x in x_ticks], rotation=45)
                    ## remove x lables
                    plt.xlabel("")
                    y_ticks = np.linspace(0, 1, 2)  
                    ax = plt.gca()
                    ax.set_yticks(y_ticks)
                    
                    pdf.savefig(bbox_inches='tight')
                    plt.close()
            #         break
            # break

def manual_main():
    out_dir = "/home/shuaiw/borg/paper/circos/borg_rev/"

    # sample = "infant_2"
    # contig = "infant_2_60_C"
    # out_dir = "/home/shuaiw/borg/paper/circos/inversion/"
    # motif_list = [["ATGCAT", 5], ["CAANNNNNNRTGA", 3], ["CAYNNNNNNTAYG", 2]]

    # sample = "soil_s4_1"
    # contig = "soil_s4_1_2981_L"

    # sample = "soil_s1_2"
    # contig = "soil_s1_2_55_L"
    # out_dir = "/home/shuaiw/borg/paper/circos/borg/"
    # motif_list = [["GATATC", 4], ["GATC", 2], ["GACGGA", 6], ["AYCNNNNNNNTCGC", 1]]

    # 


    # sample = "soil_s4_1"
    # contig = "BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GATATC", 4], ["GATC", 2], ["YCTB", 2]]

    # sample = "soil_s4_2"
    # contig = "ORANGE_SR-VP_26_10_2019_2_100cm_FINAL"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["YCT", 2], ["GCC", 2]]

    # sample = "soil_1"
    # contig = "Green_Borg_SRVP_2020_100_cm_2_FINAL-FINAL"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GGNCC", 5], ["YCTB", 2]]

    # sample = "soil_s4_2"
    # contig = "Iris_Borg_SR-VP_9_9_2021_87_5B_1_2m_2_curated"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["YCT", 2], ["ACC", 1]]

    # sample = "soil_s3_2"
    # contig = "BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GATATC", 4], ["GATC", 2], ["YCTB", 2]]

    # sample = "soil_s3_1"
    # contig = "BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GATATC", 4], ["GATC", 2], ["YCTB", 2]]

    # sample = "soil_s1_2"
    # contig = "BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GATATC", 4], ["GATC", 2], ["YCTB", 2]]

    # sample = "soil_1"
    # contig = "SRVP18_trench_6_60cm_scaf_214_117_86_FINAL"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GAA", 3],['YCTB',2],['GATC',2]]

    # sample = "soil_2"
    # contig = "SRVP18_trench_6_60cm_scaf_214_117_86_FINAL"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GAA", 3],['YCTB',2],['GATC',2]]

    sample = "soil_s1_2"
    contig = "SRVP18_trench_6_60cm_scaf_214_117_86_FINAL"
    out_dir = os.path.join(out_dir, sample)
    motif_list = [["GAA", 3],['YCTB',2],['GATC',2]]

    # sample = "soil_s3_1"
    # contig = "SRVP18_trench_6_60cm_scaf_214_117_86_FINAL"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GAA", 3],['YCTB',2],['GATC',2]]

    # sample = "soil_1"
    # contig = "Amber_Borg_SR-VP_9_9_2021_87_5B_1_2m_2_complete"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GAA", 3],['YCTB',2],['AGAA',4]]

    # sample = "soil_s4_2"
    # contig = "SR-VP_07_25_2022_A1_100cm_PACBIO-HIFI_METAMDBG_723848_L"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["YCT", 2], ["GCC", 2]]

    # sample = "soil_2"
    # contig = "SR-VP_9_9_2021_34_2B_1_4m_PACBIO-HIFI_HIFIASM-META_1061_C"
    # out_dir = os.path.join(out_dir, sample)
    # motif_list = [["GATC", 2],['AAAT',3],['GAYATC',4],['SATCAG',5],["ATGGAV",5],["CACCAG",5],["TTCGAA",6]]

    ## makdir out_dir if not exists
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    print (out_dir)
    work_dir = f"/home/shuaiw/borg/paper/borg_data/borg_for4/{sample}/{sample}_methylation3/"
    # work_dir = f"/home/shuaiw/borg/paper/run3/{sample}/{sample}_methylation3/"
    
    my_ref = f"{work_dir}/contigs/{contig}.fa"
    gff2 = f"{work_dir}/gffs/{contig}.reprocess.gff"
    ipd_ratio_file = f"{work_dir}/ipd_ratio/{contig}.ipd3.csv"
    # split_gff(gff2, contig, my_ref, out_dir)
    get_all_loci(gff2, contig, my_ref, out_dir,motif_list)

def auto_main():
    out_dir = "/home/shuaiw/borg/paper/circos/borg2/"
    all_dir = "/home/shuaiw/borg/paper/run2/"

    sample = "soil_s3_1"
    sample_obj = My_sample(sample, all_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # print (out_dir)
    # sample_obj.reference_fasta = "/home/shuaiw/borg/paper/curated_genome/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.soil_1.fa"
    # sample_obj.fai = sample_obj.reference_fasta + ".fai"
    sample_obj.read_depth()
    isolation_taxa = sample_obj.read_meta_gtdb()
    ## sort the contigs by depth
    sorted_contigs = sorted(sample_obj.depth_dict.items(), key=lambda x: x[1], reverse=True)
    for contig, depth in sorted_contigs:
        # if contig not in isolation_taxa:
        #     continue
        # if contig[-1] != "C":
        #     continue
        # taxa = isolation_taxa[contig]
        # ## check bacteria
        # if re.search("bacteria", taxa.lower()) or re.search("archaea", taxa.lower()):
        #     continue
        # if sample_obj.length_dict[contig] > 1000000:
        #     continue
        # print (contig, taxa)
        if contig != "soil_s3_1_2132_C":
            continue
        ctg_obj = My_contig(sample, all_dir, contig)
        motif_df = ctg_obj.read_motif(min_frac=0, min_sites=0)
        motif_list = []
        for _, row in motif_df.iterrows():
            if row['fraction'] > 0.1 and row["nDetected"] > 20:
                motif_list.append([row['motifString'], row['centerPos']])
            # if len(motif_list) >= 3:
            #     break

        work_dir = f"{all_dir}/{sample}/{sample}_methylation3/"
        my_ref = f"{work_dir}/contigs/{contig}.fa"
        gff2 = f"{work_dir}/gffs/{contig}.reprocess.gff"
        ipd_ratio_file = f"{work_dir}/ipd_ratio/{contig}.ipd3.csv"
        # split_gff(gff2, contig, my_ref, out_dir)
        get_all_loci(gff2, contig, my_ref, out_dir,motif_list)

def auto_main_jumbo():
    all_dir = "/home/shuaiw/borg/paper/run2/"

    out_dir = "/home/shuaiw/borg/paper/circos/jumbo/"
    jumbo_file = "/home/shuaiw/mGlu/benchmark/borg/all_jumbo_contigs_summary.tsv"

    out_dir = "/home/shuaiw/borg/paper/circos/borg3/"
    jumbo_file = "/home/shuaiw/mGlu/benchmark/borg/all_borg_contigs_summary.tsv"

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    select_jumbo = []
    for line in open(jumbo_file):
        if line.startswith("borg_ref"):
            continue
        fields = line.strip().split('\t')
        borg_ref = fields[0]
        ctg_depth = float(fields[8])
        seq_name = fields[2]
        sample_name = "_".join(seq_name.split('_')[:-2])
        if ctg_depth >= 3:
            select_jumbo.append([borg_ref, sample_name, seq_name, ctg_depth])

    for jumbo_info in select_jumbo:
        borg_ref, sample, seq_name, ctg_depth = jumbo_info
        print (borg_ref, sample, seq_name, ctg_depth)
        sample_obj = My_sample(sample, all_dir)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        sample_obj.read_depth()
        isolation_taxa = sample_obj.read_meta_gtdb()
        ## sort the contigs by depth
        sorted_contigs = sorted(sample_obj.depth_dict.items(), key=lambda x: x[1], reverse=True)
        for contig, depth in sorted_contigs:

            if contig != seq_name:
                continue
            ctg_obj = My_contig(sample, all_dir, contig)
            motif_df = ctg_obj.read_motif(min_frac=0, min_sites=0)
            if motif_df is None:
                continue
            motif_list = []
            for _, row in motif_df.iterrows():
                if row['fraction'] > 0.08 and row["nDetected"] > 50:
                    motif_list.append([row['motifString'], row['centerPos']])
                # if len(motif_list) >= 3:
                #     break

            work_dir = f"{all_dir}/{sample}/{sample}_methylation3/"
            my_ref = f"{work_dir}/contigs/{contig}.fa"
            gff2 = f"{work_dir}/gffs/{contig}.reprocess.gff"
            ipd_ratio_file = f"{work_dir}/ipd_ratio/{contig}.ipd3.csv"
            # split_gff(gff2, contig, my_ref, out_dir)
            get_all_loci(gff2, contig, my_ref, out_dir,motif_list, depth=depth)


def auto_main_borg():
    out_dir = "/home/shuaiw/borg/paper/circos/mini_borg/"
    all_dir = "/home/shuaiw/borg/paper/borg_data/mini_borg/"

    ## makdir out_dir if not exists


    sample = "soil_1"
    sample_obj = My_sample(sample, all_dir)
    out_dir = os.path.join(out_dir, sample)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    sample_obj.reference_fasta = "/home/shuaiw/borg/paper/borg_data/mp_ece.contigs.fa"
    sample_obj.fai = sample_obj.reference_fasta + ".fai"
    sample_obj.read_depth()

    sorted_contigs = sorted(sample_obj.depth_dict.items(), key=lambda x: x[1], reverse=True)
    for contig, depth in sorted_contigs:

        if contig not in ["Saturn_mini-Borg_41kb_SR-VP_26_10_2019_2_90cm_33_22",\
                          "Uranus_mini-Borg_55kb_SR-VP_9_9_2021_81_5A_0_75m_35_16_complete",\
                          "SR-VP_26_10_2019_1_90cm_scaffold_8229",\
                          "SR-VP_26_10_2019_2_90cm_scaffold_29557",\
                            "Saturn_mini-Borg_75kb_SRVP18_trench_6_60cm_34_20_near_complete",\
                                "SR-VP_26_10_2019_2_90cm_scaffold_28500",\
                                    "SR-VP_0-2_scaffold_141_3651159",\
                                        "SR-VP_26_10_2020_1_100CM_scaffold_3695492",\
                                            "SRVP18_hole-7m-from-trench_1_20cm_scaffold_8051",\
                                                "SRVP18_hole-7m-from-trench_1_20cm_scaffold_3869",\
                                                    "SRVP18_hole-7m-from-trench_1_20cm_scaffold_1996",\
                                                        "SR-VP_26_10_2019_2_100cm_scaffold_2211", \
                                                            "SR-VP_26_10_2020_1_100CM_scaffold_6826853",\
                                                                "Uranus_mini-Borg_69kb_SR-VP_26_10_2020_2_100CM_33_33",\
                                                                    "Saturn_mini-Borg_94kb_SR-VP_9_9_2021_87_5B_1_2m_33_10_near_complete",\
                                                                        "Venus_mini-Borg_106kb_SRVP18_hole-7m-from-trench_1_80cm_34_23_complete"]:
            continue

        ctg_obj = My_contig(sample, all_dir, contig)
        motif_df = ctg_obj.read_motif(min_frac=0, min_sites=0)
        if motif_df is None:
            continue
        motif_list = []
        for _, row in motif_df.iterrows():
            if row['fraction'] > 0.1 and row["nDetected"] > 20:
                motif_list.append([row['motifString'], row['centerPos']])
        print (motif_list)
        work_dir = f"{all_dir}/{sample}/{sample}_methylation3/"
        my_ref = f"{work_dir}/contigs/{contig}.fa"
        gff2 = f"{work_dir}/gffs/{contig}.reprocess.gff"

        get_all_loci(gff2, contig, my_ref, out_dir,motif_list)

if __name__ == "__main__":
    score_cutoff = 30
    # auto_main()
    # auto_main_jumbo()
    manual_main()
    # auto_main_borg()


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



    




