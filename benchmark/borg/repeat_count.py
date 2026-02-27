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
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import get_unique_motifs, My_sample, Isolation_sample, My_contig, My_cluster, classify_taxa, get_ctg_taxa,get_detail_taxa_name
from scipy.stats import fisher_exact
# ...existing code...
from scipy.stats import fisher_exact

def read_ref(ref):
    REF = {}
    for record in SeqIO.parse(ref, "fasta"):
    #     seq_dict[record.id] = record.seq
    # return seq_dict
        REF[record.id] = record.seq
        # return str(record.seq), record.id
    return REF

def get_motif_sites(REF, motif_new, exact_pos, modified_loci, max_len = 1000000000):
    motif_len = len(motif_new)
    rev_exact_pos = motif_len - exact_pos + 1
    motif_sites = {}
    motif_loci_num = 0
    motif_modify_num = 0
    for_loci_num = 0
    rev_loci_num = 0
    for_modified_num = 0
    rev_modified_num = 0

    record_modified_sites = {}
    record_motif_sites = defaultdict(lambda: {"+":[], "-":[]})
    record_mod_motif_sites = defaultdict(lambda: {"+":[], "-":[]})

    for r, contig in REF.items():
        contig = contig[:max_len]
        for site in nt_search(str(contig), motif_new)[1:]:
            tag = r + ":" + str(site+exact_pos) + "+"
            record_modified_sites[tag] = motif_new
            record_motif_sites[r]["+"].append(site+exact_pos)
            motif_loci_num += 1
            for_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                for_modified_num += 1
                record_mod_motif_sites[r]["+"].append(site+exact_pos)

        for site in nt_search(str(contig), Seq(motif_new).reverse_complement())[
            1:
        ]:
            tag = r + ":" + str(site+rev_exact_pos) + "-"
            record_modified_sites[tag] = motif_new
            record_motif_sites[r]["-"].append(site+rev_exact_pos)
            motif_loci_num += 1
            rev_loci_num += 1
            if tag in modified_loci:
                motif_modify_num += 1
                rev_modified_num += 1
                record_mod_motif_sites[r]["-"].append(site+rev_exact_pos)
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

    print (for_loci_num, for_modified_num,for_ratio,\
            rev_loci_num, rev_modified_num, rev_ratio,\
            motif_loci_num, motif_modify_num, ratio, proportion_all_modified)
    return [for_loci_num, for_modified_num,for_ratio,\
            rev_loci_num, rev_modified_num, rev_ratio,\
            motif_loci_num, motif_modify_num, ratio, proportion_all_modified], record_modified_sites, record_motif_sites, record_mod_motif_sites

def get_modified_ratio(gff, score_cutoff = 30):
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

def read_repeat_bed(repeat_bed):
    """Read repeat regions from BED file"""
    repeat_regions = defaultdict(list)
    with open(repeat_bed, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip comments
                continue
            fields = line.strip().split("\t")
            seq_id = fields[0]
            start = int(fields[1])  # BED format: 0-based start
            end = int(fields[2])    # BED format: exclusive end
            repeat_regions[seq_id].append((start, end))
    return repeat_regions

def read_orf_gff(gff_file):
    """Read ORF regions from GFF file"""
    orf_regions = defaultdict(list)
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip comments
                continue
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            # Only process ORF features
            if feature_type == 'ORF':
                seq_id = fields[0]
                start = int(fields[3])  # GFF format: 1-based start
                end = int(fields[4])    # GFF format: inclusive end
                orf_regions[seq_id].append((start, end))
    return orf_regions

def is_in_repeat(pos, repeat_regions):
    """Check if a position is within any repeat region"""
    for start, end in repeat_regions:
        if start <= pos <= end:
            return True
    return False

def is_in_region(pos, regions):
    """Check if a position is within any region (generic function)"""
    for start, end in regions:
        if start <= pos <= end:
            return True
    return False

def analyze_motif_repeat_enrichment(REF, record_motif_sites, repeat_regions):
    """Count motif sites in and out of repeat regions"""
    in_repeat = 0
    out_repeat = 0
    
    for contig_id, motif_sites in record_motif_sites.items():
        if contig_id not in repeat_regions:
            # No repeats in this contig, all motifs are outside
            out_repeat += len(motif_sites['+']) + len(motif_sites['-'])
            continue
            
        repeats = repeat_regions[contig_id]
        
        # Check forward strand
        for pos in motif_sites['+']:
            if is_in_repeat(pos, repeats):
                in_repeat += 1
            else:
                out_repeat += 1
        
        # Check reverse strand
        for pos in motif_sites['-']:
            if is_in_repeat(pos, repeats):
                in_repeat += 1
            else:
                out_repeat += 1
    
    # Calculate total repeat and non-repeat bases
    repeat_bases = 0
    for contig_id, repeats in repeat_regions.items():
        for start, end in repeats:
            repeat_bases += (end - start + 1)
    
    contig_total_bases = sum(len(seq) for seq in REF.values())
    non_repeat_bases = contig_total_bases - repeat_bases
    
    # Fisher's exact test
    # Contingency table: [[in_repeat, out_repeat], [repeat_bases - in_repeat, non_repeat_bases - out_repeat]]
    contingency_table = [[in_repeat, out_repeat], 
                         [repeat_bases - in_repeat, non_repeat_bases - out_repeat]]
    
    oddsratio, pvalue = fisher_exact(contingency_table)
    
    print(f"\n{'='*60}")
    print(f"MOTIF ENRICHMENT IN REPEAT REGIONS")
    print(f"{'='*60}")
    print(f"Motif sites in repeat regions: {in_repeat}")
    print(f"Motif sites outside repeat regions: {out_repeat}")
    print(f"Total repeat bases: {repeat_bases:,}")
    print(f"Total non-repeat bases: {non_repeat_bases:,}")
    print(f"Total bases: {contig_total_bases:,}")
    print(f"Repeat coverage: {repeat_bases/contig_total_bases*100:.2f}%")
    print(f"\nFisher's exact test:")
    print(f"Odds ratio: {oddsratio:.4f}")
    print(f"P-value: {pvalue:.4e}")
    
    if pvalue < 0.05:
        if oddsratio > 1:
            print("Result: Motif sites are significantly ENRICHED in repeat regions ✓")
        else:
            print("Result: Motif sites are significantly DEPLETED in repeat regions ✗")
    else:
        print("Result: No significant enrichment or depletion")
    print(f"{'='*60}\n")
    
    return in_repeat, out_repeat, repeat_bases, non_repeat_bases, oddsratio, pvalue

def analyze_motif_orf_enrichment(REF, record_motif_sites, orf_regions):
    """Count motif sites in and out of ORF regions"""
    in_orf = 0
    out_orf = 0
    
    for contig_id, motif_sites in record_motif_sites.items():
        if contig_id not in orf_regions:
            # No ORFs in this contig, all motifs are outside
            out_orf += len(motif_sites['+']) + len(motif_sites['-'])
            continue
            
        orfs = orf_regions[contig_id]
        
        # Check forward strand
        for pos in motif_sites['+']:
            if is_in_region(pos, orfs):
                in_orf += 1
            else:
                out_orf += 1
        
        # Check reverse strand
        for pos in motif_sites['-']:
            if is_in_region(pos, orfs):
                in_orf += 1
            else:
                out_orf += 1
    
    # Calculate total ORF and non-ORF bases
    orf_bases = 0
    for contig_id, orfs in orf_regions.items():
        for start, end in orfs:
            orf_bases += (end - start + 1)
    
    contig_total_bases = sum(len(seq) for seq in REF.values())
    non_orf_bases = contig_total_bases - orf_bases
    
    # Fisher's exact test
    # Contingency table: [[in_orf, out_orf], [orf_bases - in_orf, non_orf_bases - out_orf]]
    contingency_table = [[in_orf, out_orf], 
                         [orf_bases - in_orf, non_orf_bases - out_orf]]
    
    oddsratio, pvalue = fisher_exact(contingency_table)
    
    print(f"\n{'='*60}")
    print(f"MOTIF ENRICHMENT IN ORF REGIONS")
    print(f"{'='*60}")
    print(f"Motif sites in ORF regions: {in_orf}")
    print(f"Motif sites outside ORF regions: {out_orf}")
    print(f"Total ORF bases: {orf_bases:,}")
    print(f"Total non-ORF bases: {non_orf_bases:,}")
    print(f"Total bases: {contig_total_bases:,}")
    print(f"ORF coverage: {orf_bases/contig_total_bases*100:.2f}%")
    print(f"\nFisher's exact test:")
    print(f"Odds ratio: {oddsratio:.4f}")
    print(f"P-value: {pvalue:.4e}")
    
    if pvalue < 0.05:
        if oddsratio > 1:
            print("Result: Motif sites are significantly ENRICHED in ORF regions ✓")
        else:
            print("Result: Motif sites are significantly DEPLETED in ORF regions ✗")
    else:
        print("Result: No significant enrichment or depletion")
    print(f"{'='*60}\n")
    
    return in_orf, out_orf, orf_bases, non_orf_bases, oddsratio, pvalue




ORF_anno = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.gff"
mod_gff = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR_score30.gff"
borg_ref = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.fasta"
repeat_bed = "/home/shuaiw/borg/paper/borg_data/batch_export2/repeat.bed"
motif_new = "YCTK"
exact_pos = 2

modified_loci = get_modified_ratio(mod_gff, score_cutoff = 30)
REF = read_ref(borg_ref)
contig_len = sum(len(seq) for seq in REF.values())
x, y, record_motif_sites, record_mod_motif_sites = get_motif_sites(REF, motif_new, exact_pos, modified_loci)

# Analyze motif enrichment in repeat regions
repeat_regions = read_repeat_bed(repeat_bed)
analyze_motif_repeat_enrichment(REF, record_motif_sites, repeat_regions)
analyze_motif_repeat_enrichment(REF, record_mod_motif_sites, repeat_regions)

# Analyze motif enrichment in ORF regions
orf_regions = read_orf_gff(ORF_anno)
analyze_motif_orf_enrichment(REF, record_motif_sites, orf_regions)
analyze_motif_orf_enrichment(REF, record_mod_motif_sites, orf_regions)