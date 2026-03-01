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
            if feature_type == 'ORF' or feature_type == 'CDS':
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

def compare_more(REF, record_motif_sites, orf_regions, repeat_regions):
    """Compare enrichment across three region types: ORF, Repeat, and Intergenic (non-ORF & non-repeat)"""
    
    in_orf_only = 0
    in_repeat_only = 0  
    in_intergenic = 0  # neither ORF nor repeat
    in_both = 0  # in both ORF and repeat (overlap)
    
    for contig_id, motif_sites in record_motif_sites.items():
        orfs = orf_regions.get(contig_id, [])
        repeats = repeat_regions.get(contig_id, [])
        
        # Check all motif positions (both strands)
        all_positions = motif_sites['+'] + motif_sites['-']
        
        for pos in all_positions:
            in_orf = is_in_region(pos, orfs)
            in_repeat = is_in_region(pos, repeats)
            
            if in_orf and in_repeat:
                in_both += 1
            elif in_orf:
                in_orf_only += 1
            elif in_repeat:
                in_repeat_only += 1
            else:
                in_intergenic += 1
    
    # Calculate base counts for each region type
    orf_bases = 0
    repeat_bases = 0
    overlap_bases = 0
    
    for contig_id in REF.keys():
        orfs = orf_regions.get(contig_id, [])
        repeats = repeat_regions.get(contig_id, [])
        
        # Count ORF bases
        for start, end in orfs:
            orf_bases += (end - start + 1)
        
        # Count repeat bases
        for start, end in repeats:
            repeat_bases += (end - start + 1)
        
        # Count overlap between ORF and repeat
        for orf_start, orf_end in orfs:
            for rep_start, rep_end in repeats:
                # Calculate overlap
                overlap_start = max(orf_start, rep_start)
                overlap_end = min(orf_end, rep_end)
                if overlap_start <= overlap_end:
                    overlap_bases += (overlap_end - overlap_start + 1)
    
    contig_total_bases = sum(len(seq) for seq in REF.values())
    orf_only_bases = orf_bases - overlap_bases
    repeat_only_bases = repeat_bases - overlap_bases
    intergenic_bases = contig_total_bases - orf_bases - repeat_bases + overlap_bases
    
    # Fisher's exact test: ORF-only vs Intergenic
    contingency_orf_vs_intergenic = [[in_orf_only, in_intergenic], 
                                      [orf_only_bases - in_orf_only, intergenic_bases - in_intergenic]]
    odds_orf_intergenic, pval_orf_intergenic = fisher_exact(contingency_orf_vs_intergenic)
    
    # Fisher's exact test: Repeat-only vs Intergenic
    contingency_rep_vs_intergenic = [[in_repeat_only, in_intergenic],
                                      [repeat_only_bases - in_repeat_only, intergenic_bases - in_intergenic]]
    odds_rep_intergenic, pval_rep_intergenic = fisher_exact(contingency_rep_vs_intergenic)
    
    # Fisher's exact test: ORF-only vs Repeat-only
    contingency_orf_vs_rep = [[in_orf_only, in_repeat_only],
                               [orf_only_bases - in_orf_only, repeat_only_bases - in_repeat_only]]
    odds_orf_rep, pval_orf_rep = fisher_exact(contingency_orf_vs_rep)
    
    print(f"\n{'='*60}")
    print(f"COMPARATIVE ENRICHMENT ANALYSIS")
    print(f"{'='*60}")
    print(f"\nMotif Distribution:")
    print(f"  ORF-only regions: {in_orf_only}")
    print(f"  Repeat-only regions: {in_repeat_only}")
    print(f"  Intergenic regions (non-ORF & non-repeat): {in_intergenic}")
    print(f"  Overlap (ORF & Repeat): {in_both}")
    print(f"  Total motifs: {in_orf_only + in_repeat_only + in_intergenic + in_both}")
    
    print(f"\nBase Distribution:")
    print(f"  ORF-only bases: {orf_only_bases:,} ({orf_only_bases/contig_total_bases*100:.2f}%)")
    print(f"  Repeat-only bases: {repeat_only_bases:,} ({repeat_only_bases/contig_total_bases*100:.2f}%)")
    print(f"  Intergenic bases: {intergenic_bases:,} ({intergenic_bases/contig_total_bases*100:.2f}%)")
    print(f"  Overlap bases: {overlap_bases:,} ({overlap_bases/contig_total_bases*100:.2f}%)")
    
    print(f"\nMotif Density (motifs per kb):")
    print(f"  ORF-only: {in_orf_only/(orf_only_bases/1000):.2f}")
    print(f"  Repeat-only: {in_repeat_only/(repeat_only_bases/1000) if repeat_only_bases > 0 else 0:.2f}")
    print(f"  Intergenic: {in_intergenic/(intergenic_bases/1000):.2f}")
    
    print(f"\nStatistical Comparisons (Fisher's Exact Test):")
    print(f"\n  ORF-only vs Intergenic:")
    print(f"    Odds ratio: {odds_orf_intergenic:.4f}")
    print(f"    P-value: {pval_orf_intergenic:.4e}")
    if pval_orf_intergenic < 0.05:
        result = "ENRICHED" if odds_orf_intergenic > 1 else "DEPLETED"
        print(f"    Result: ORFs are significantly {result} vs intergenic regions")

def compare_TR(REF, orf_regions, repeat_regions):
    """Compare the enrichment of tandem repeats in ORF regions vs non-ORF regions"""
    
    repeat_in_orf = 0
    repeat_in_non_orf = 0
    
    for contig_id in REF.keys():
        orfs = orf_regions.get(contig_id, [])
        repeats = repeat_regions.get(contig_id, [])
        
        # For each repeat region, check if it overlaps with ORFs
        for rep_start, rep_end in repeats:
            repeat_bases = rep_end - rep_start + 1
            overlap_with_orf = 0
            
            # Check overlap with each ORF
            for orf_start, orf_end in orfs:
                overlap_start = max(rep_start, orf_start)
                overlap_end = min(rep_end, orf_end)
                if overlap_start <= overlap_end:
                    overlap_with_orf += (overlap_end - overlap_start + 1)
            
            # Count bases in ORF vs non-ORF
            repeat_in_orf += overlap_with_orf
            repeat_in_non_orf += (repeat_bases - overlap_with_orf)
    
    # Calculate total ORF and non-ORF bases
    orf_bases = 0
    for contig_id, orfs in orf_regions.items():
        for start, end in orfs:
            orf_bases += (end - start + 1)
    
    contig_total_bases = sum(len(seq) for seq in REF.values())
    non_orf_bases = contig_total_bases - orf_bases
    
    # Calculate non-repeat bases in each region
    non_repeat_in_orf = orf_bases - repeat_in_orf
    non_repeat_in_non_orf = non_orf_bases - repeat_in_non_orf
    
    # Fisher's exact test
    contingency_table = [[repeat_in_orf, repeat_in_non_orf],
                         [non_repeat_in_orf, non_repeat_in_non_orf]]
    
    oddsratio, pvalue = fisher_exact(contingency_table)
    
    # Calculate densities
    repeat_density_orf = (repeat_in_orf / orf_bases * 100) if orf_bases > 0 else 0
    repeat_density_non_orf = (repeat_in_non_orf / non_orf_bases * 100) if non_orf_bases > 0 else 0
    
    print(f"\n{'='*60}")
    print(f"TANDEM REPEAT ENRICHMENT IN ORF vs NON-ORF REGIONS")
    print(f"{'='*60}")
    print(f"\nRepeat Distribution:")
    print(f"  Repeat bases in ORF regions: {repeat_in_orf:,}")
    print(f"  Repeat bases in non-ORF regions: {repeat_in_non_orf:,}")
    print(f"  Total repeat bases: {repeat_in_orf + repeat_in_non_orf:,}")
    
    print(f"\nRegion Sizes:")
    print(f"  Total ORF bases: {orf_bases:,} ({orf_bases/contig_total_bases*100:.2f}%)")
    print(f"  Total non-ORF bases: {non_orf_bases:,} ({non_orf_bases/contig_total_bases*100:.2f}%)")
    print(f"  Total bases: {contig_total_bases:,}")
    
    print(f"\nRepeat Density:")
    print(f"  In ORF regions: {repeat_density_orf:.4f}% of ORF bases are repeats")
    print(f"  In non-ORF regions: {repeat_density_non_orf:.4f}% of non-ORF bases are repeats")
    print(f"  Fold enrichment: {repeat_density_orf/repeat_density_non_orf if repeat_density_non_orf > 0 else 'N/A':.2f}x")
    
    print(f"\nFisher's exact test:")
    print(f"Odds ratio: {oddsratio:.4f}")
    print(f"P-value: {pvalue:.4e}")
    
    if pvalue < 0.05:
        if oddsratio > 1:
            print("Result: Tandem repeats are significantly ENRICHED in ORF regions ✓")
        else:
            print("Result: Tandem repeats are significantly DEPLETED in ORF regions ✗")
    else:
        print("Result: No significant enrichment or depletion")
    print(f"{'='*60}\n")
    
    return repeat_in_orf, repeat_in_non_orf, orf_bases, non_orf_bases, oddsratio, pvalue

def compare_motif_in_TR(REF, record_motif_sites, orf_regions, repeat_regions, Modflag = "no"):
    """Compare motif occurrence in TRs within ORF regions vs TRs in intergenic regions"""
    
    # Count motifs in TRs that are in ORF regions vs intergenic regions
    motif_in_TR_in_orf = 0
    motif_in_TR_in_intergenic = 0
    
    # Count TR bases in each region type
    TR_bases_in_orf = 0
    TR_bases_in_intergenic = 0
    
    for contig_id, motif_sites in record_motif_sites.items():
        orfs = orf_regions.get(contig_id, [])
        repeats = repeat_regions.get(contig_id, [])
        # print (Modflag, contig_id)
        # Check all motif positions (both strands)
        all_positions = motif_sites['+'] + motif_sites['-']
        
        for pos in all_positions:
            in_repeat = is_in_region(pos, repeats)
            if in_repeat:
                # This motif is in a TR, now check if the TR is in ORF or intergenic
                in_orf = is_in_region(pos, orfs)
                if in_orf:
                    motif_in_TR_in_orf += 1
                else:
                    motif_in_TR_in_intergenic += 1
    
    # Calculate TR bases in ORF vs intergenic regions
    for contig_id in REF.keys():
        orfs = orf_regions.get(contig_id, [])
        repeats = repeat_regions.get(contig_id, [])
        
        # For each repeat region, determine how much overlaps with ORF
        for rep_start, rep_end in repeats:
            repeat_bases = rep_end - rep_start + 1
            overlap_with_orf = 0
            
            # Check overlap with each ORF
            for orf_start, orf_end in orfs:
                overlap_start = max(rep_start, orf_start)
                overlap_end = min(rep_end, orf_end)
                if overlap_start <= overlap_end:
                    overlap_with_orf += (overlap_end - overlap_start + 1)
            
            # Count TR bases in ORF vs intergenic
            TR_bases_in_orf += overlap_with_orf
            TR_bases_in_intergenic += (repeat_bases - overlap_with_orf)
    
    # Calculate non-motif TR bases
    non_motif_TR_in_orf = TR_bases_in_orf - motif_in_TR_in_orf
    non_motif_TR_in_intergenic = TR_bases_in_intergenic - motif_in_TR_in_intergenic
    
    # Fisher's exact test
    contingency_table = [[motif_in_TR_in_orf, motif_in_TR_in_intergenic],
                         [non_motif_TR_in_orf, non_motif_TR_in_intergenic]]
    
    oddsratio, pvalue = fisher_exact(contingency_table)
    
    # Calculate densities
    density_TR_in_orf = (motif_in_TR_in_orf / TR_bases_in_orf * 1000) if TR_bases_in_orf > 0 else 0
    density_TR_in_intergenic = (motif_in_TR_in_intergenic / TR_bases_in_intergenic * 1000) if TR_bases_in_intergenic > 0 else 0
    
    print(f"\n{'='*60}")
    print(f"MOTIF ENRICHMENT IN TANDEM REPEATS:")
    print(f"ORF-LOCATED TRs vs INTERGENIC-LOCATED TRs")
    print(f"{'='*60}")
    print(f"\nMotif Distribution in TRs:")
    print(f"  Motifs in TRs within ORF regions: {motif_in_TR_in_orf}")
    print(f"  Motifs in TRs within intergenic regions: {motif_in_TR_in_intergenic}")
    print(f"  Total motifs in TRs: {motif_in_TR_in_orf + motif_in_TR_in_intergenic}")
    
    print(f"\nTR Base Distribution:")
    print(f"  TR bases in ORF regions: {TR_bases_in_orf:,}")
    print(f"  TR bases in intergenic regions: {TR_bases_in_intergenic:,}")
    print(f"  Total TR bases: {TR_bases_in_orf + TR_bases_in_intergenic:,}")
    
    print(f"\nMotif Density in TRs (motifs per kb):")
    print(f"  TRs in ORF regions: {density_TR_in_orf:.2f}")
    print(f"  TRs in intergenic regions: {density_TR_in_intergenic:.2f}")
    if TR_bases_in_intergenic > 0 and density_TR_in_intergenic > 0:
        print(f"  Fold difference: {density_TR_in_orf/density_TR_in_intergenic:.2f}x")
    
    print(f"\nFisher's exact test:")
    print(f"Odds ratio: {oddsratio:.4f}")
    print(f"P-value: {pvalue:.4e}")
    
    if pvalue < 0.05:
        if oddsratio > 1:
            print("Result: Motifs are significantly MORE abundant in ORF-located TRs ✓")
        else:
            print("Result: Motifs are significantly LESS abundant in ORF-located TRs ✗")
    else:
        print("Result: No significant difference between ORF-located and intergenic TRs")
    print(f"{'='*60}\n")
    
    return motif_in_TR_in_orf, motif_in_TR_in_intergenic, TR_bases_in_orf, TR_bases_in_intergenic, oddsratio, pvalue

if __name__ == "__main__":

    # ORF_anno = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.gff"
    # mod_gff = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR_score30.gff"
    # borg_ref = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.fasta"
    # repeat_bed = "/home/shuaiw/borg/paper/borg_data/batch_export2/repeat.bed"

    folder = "/home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/"

    # ORF_anno = f"{folder}/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.prodigal.gff"
    # mod_gff = f"/home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_115/gffs/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.all.reprocess.gff"
    # borg_ref = f"{folder}/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.contigs.fa"
    # repeat_bed = f"{folder}/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.repeat.bed"

    # ORF_anno = f"{folder}/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.prodigal.gff"
    # mod_gff = f"/home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_90/gffs/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.gff"
    # borg_ref = f"{folder}/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.fa"
    # repeat_bed = f"{folder}/SR-VP_07_25_2022_A1_90cm_PACBIO-HIFI_METAMDBG_641677_L.repeat.bed"

    # ORF_anno = f"{folder}/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_almost_complete_Black_borg_32_00.prodigal.gff"
    # mod_gff = f"/home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_80/gffs/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_METAMDBG_428814_L.gff"
    # borg_ref = f"{folder}/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_almost_complete_Black_borg_32_00.contigs.fa"
    # repeat_bed = f"{folder}/SR-VP_07_25_2022_A1_80cm_PACBIO-HIFI_almost_complete_Black_borg_32_00.repeat.bed"

    ORF_anno = f"{folder}/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.prodigal.gff"
    mod_gff = f"/home/shuaiw/borg/paper/borg_data/batch_export2/new_run/soil_60/gffs/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.all.reprocess.gff"
    borg_ref = f"{folder}/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.contigs.fa"
    repeat_bed = f"{folder}/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_Black_Borg_32_00.repeat.bed"



    host_motifs_file = "/home/shuaiw/borg/paper/gg_run3/soil_60/soil_60_methylation4/motifs/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_METAMDBG_551173_L.motifs.csv"
    
    
    output_dir="../../tmp/figures/borg_fig"
    motif_new = "YCTK"
    # motif_new = "GATC"
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

    # Compare enrichment across ORF, Repeat, and Intergenic regions
    print("\n### Analysis for ALL motif sites ###")
    compare_more(REF, record_motif_sites, orf_regions, repeat_regions)

    print("\n### Analysis for MODIFIED motif sites only ###")
    compare_more(REF, record_mod_motif_sites, orf_regions, repeat_regions)

    # Compare tandem repeat enrichment in ORF vs non-ORF regions
    compare_TR(REF, orf_regions, repeat_regions)

    # Compare motif occurrence in TRs within ORFs vs TRs in intergenic regions
    print("\n### Analysis for ALL motif sites in TRs ###")
    compare_motif_in_TR(REF, record_motif_sites, orf_regions, repeat_regions)

    print("\n### Analysis for MODIFIED motif sites in TRs ###")
    compare_motif_in_TR(REF, record_mod_motif_sites, orf_regions, repeat_regions, Modflag="yes")


