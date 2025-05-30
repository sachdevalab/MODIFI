from Bio.SeqUtils import nt_search
from Bio import SeqIO
from Bio.Seq import Seq
import xml.etree.ElementTree as ET
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import sys
from scipy.stats import fisher_exact
import os

def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        # return str(record.seq), record.id
    return REF

def if_region_in_gene_regulatory_region(pos, Gene_regulatory_regions, strand):
    for start, end, cds_strand in Gene_regulatory_regions:
        if cds_strand == strand and start <= pos <= end:
            return True
    return False

def get_motif_sites(REF, motif_new, exact_pos, modified_loci, Gene_regulatory_regions, background_regions):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos + 1

    gene_regulatory_region_num = 0
    methylated_gene_regulatory_region_num = 0
    normal_region_num = 0
    methlated_normal_region_num = 0

    control_region_num = 0
    methylated_control_num = 0

    record_modified_sites = {}

    for r, contig in REF.items():
        for site in nt_search(str(contig), motif_new)[1:]:

            tag = r + ":" + str(site+exact_pos) + "+"
            ## check if the site is in the gene regulatory regions
            if if_region_in_gene_regulatory_region(site + exact_pos, Gene_regulatory_regions, "+"):
                gene_regu_flag = True
                gene_regulatory_region_num += 1
            else:
                gene_regu_flag = False
                normal_region_num += 1

            if tag in modified_loci:
                if gene_regu_flag:
                    methylated_gene_regulatory_region_num += 1
                else:
                    methlated_normal_region_num += 1

            if if_region_in_gene_regulatory_region(site + exact_pos, background_regions, "+"):
                control_region_num += 1
                if tag in modified_loci:
                    methylated_control_num += 1


        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            tag = r + ":" + str(site+rev_exact_pos) + "-"
            record_modified_sites[tag] = motif_new


            if if_region_in_gene_regulatory_region(site + rev_exact_pos, Gene_regulatory_regions, "-"):
                gene_regu_flag = True
                gene_regulatory_region_num += 1
            else:
                gene_regu_flag = False
                normal_region_num += 1
            if tag in modified_loci:
                if gene_regu_flag:
                    methylated_gene_regulatory_region_num += 1
                else:
                    methlated_normal_region_num += 1

            if if_region_in_gene_regulatory_region(site + exact_pos, background_regions, "-"):
                control_region_num += 1
                if tag in modified_loci:
                    methylated_control_num += 1

    # print ("gene_regulatory_region_num", gene_regulatory_region_num, "methylated_gene_regulatory_region_num", methylated_gene_regulatory_region_num, "methylated ratio", methylated_gene_regulatory_region_num/gene_regulatory_region_num)
    # print ("normal_region_num", normal_region_num, "methlated_normal_region_num", methlated_normal_region_num, "methylated ratio", methlated_normal_region_num/normal_region_num)
    ## fisher test for the two proportions
    # a, b, c, d = methylated_gene_regulatory_region_num, gene_regulatory_region_num - methylated_gene_regulatory_region_num, methlated_normal_region_num, normal_region_num - methlated_normal_region_num
    a, b, c, d = methylated_gene_regulatory_region_num, gene_regulatory_region_num - methylated_gene_regulatory_region_num, methylated_control_num, control_region_num - methylated_control_num
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
    # print ("odds ratio", odds_ratio, "p-value", p_value)
    return odds_ratio, p_value, a, b, c, d

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
    # print ("no. of modified loci", len(modified_loci))
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


def collect_regulation_region(gff, genome):
    Gene_regulatory_regions  = []
    for line in open(gff, "r"):
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        if len(line) < 9:
            continue
        if line[0] != genome:
            continue
        if line[2] == "CDS":
            # print (line)
            start = int(line[3])
            end = int(line[4])
            strand = line[6]
            #Gene regulatory regions were defined as 100bp upstream of the CDS start to 50bp downstream of the CDS start. 
            if strand == "+":
                reg_start = max(1, start - 100)
                reg_end = start + 50
                Gene_regulatory_regions.append((reg_start, reg_end, "+"))
            elif strand == "-":
                reg_start = max(1, end - 50)
                reg_end = end + 100
                Gene_regulatory_regions.append((reg_start, reg_end, "-"))
    # print (len(Gene_regulatory_regions), "gene regulatory regions collected")
    return Gene_regulatory_regions

def get_background(Gene_regulatory_regions):
    max_len = Gene_regulatory_regions[-1][1]
    ## randomly select 1000 regions from the genome, 500 in + strand and 500 in - strand
    import random
    region_length = 150
    background_regions = []
    for _ in range(500):
        start = random.randint(1, max_len - region_length)
        end = start + region_length
        background_regions.append((start, end, "+"))
    for _ in range(500):
        start = random.randint(1, max_len - region_length)
        end = start + region_length
        background_regions.append((start, end, "-"))
    return background_regions


if __name__ == "__main__":
    score_cutoff = 30
    genome_gff = "/home/shuaiw/borg/pengfan/contigs/protein/Final_Genomes_qc_rmcirc_prodigal_complete.gff"
    genome = "RuReacBro_20230708_10_40h_50ppm_r1_scaffold_3"

    data = []
    for file in os.listdir("/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/motifs/"):
        if file.endswith(".motifs.csv"):
            ## check if the file has more than one line
            with open(f"/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/motifs/{file}", "r") as f:
                lines = f.readlines()
                if len(lines) < 2:
                    # print(f"File {file} has less than two lines, skipping.")
                    continue
            genome = file.split(".")[0]

            my_ref = f"/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/contigs/{genome}.fa"
            gff = f"/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/gffs/{genome}.gff"
            motif_file = f"/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/motifs/{genome}.motifs.csv"

            # motif_new = "CAGNNNNNNTRG"
            # exact_pos = 2


            Gene_regulatory_regions = collect_regulation_region(genome_gff, genome)
            if len(Gene_regulatory_regions) == 0:
                continue
            background_regions = get_background(Gene_regulatory_regions)

            REF = read_ref(my_ref)
            # print (REF)
            modified_loci = get_modified_ratio(gff)

            motifs = pd.read_csv(motif_file)

            for index, motif in motifs.iterrows():
                motif_new = motif["motifString"]
                exact_pos = motif["centerPos"]
                odds_ratio, p_value, a, b, c, d = get_motif_sites(REF, motif_new, exact_pos, modified_loci, Gene_regulatory_regions, background_regions)
                if p_value <= 0.05:
                    print(f"genome {genome} , Motif: {motif_new}, Odds Ratio: {odds_ratio}, P-value: {p_value}", a, b, c, d)
                    data.append([genome, motif_new, exact_pos, odds_ratio, p_value, a, b, c, d])
    df = pd.DataFrame(data, columns=["genome", "motif", "exact_pos", "odds_ratio", "p_value", "a", "b", "c", "d"])
    df.to_csv("/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin2/regulatory_motif_enrichment.csv", index=False)



# motifString,centerPos,modificationType,fraction,nDetected,nGenome,groupTag,partnerMotifString,meanScore,meanIpdRatio,meanCoverage,objectiveScore
# CYAABTC,4,modified_base,0.7329193,118,161,CYAABTC,,50.966103,1.534779,80.60169,4546.885
# CAANNNNNGTTG,3,modified_base,0.7446808,35,47,CAANNNNNGTTG,,45.885715,1.4865172,76.71429,1231.739
# CAGNNNNNNTRG,2,modified_base,0.3181818,70,220,CAGNNNNNNTRG,,45.75714,1.4647014,76.6,1142.7856
