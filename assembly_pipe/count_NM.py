import pysam
import os
from multiprocessing import Pool
import sys
import argparse
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re

def parse_cigar_operations(cigar_tuples):
    """
    Parse CIGAR operations to count matches, mismatches, insertions, and deletions
    
    CIGAR operations:
    0: M (match/mismatch)
    1: I (insertion)
    2: D (deletion)
    4: S (soft clipping)
    5: H (hard clipping)
    7: = (sequence match)
    8: X (sequence mismatch)
    """
    matches = 0
    mismatches = 0
    insertions = 0
    deletions = 0
    soft_clip = 0
    hard_clip = 0
    i = 0
    cigar_num =len(cigar_tuples)
    for op, length in cigar_tuples:
        if op == 0:  # M (match/mismatch) - need NM tag to distinguish
            matches += length  # We'll adjust this with NM tag
        elif op == 1:  # I (insertion)
            insertions += length
        elif op == 2:  # D (deletion)
            deletions += length
        elif op == 7:  # = (sequence match)
            matches += length
        elif op == 8:  # X (sequence mismatch)
            mismatches += length
        elif op == 4:  # S (soft clipping)
            if i != 0 and i != cigar_num - 1:   ### ignore clip at ends
                soft_clip += length
        elif op == 5:  # H (hard clipping)
            if i != 0 and i != cigar_num - 1:   ### ignore clip at ends
                hard_clip += length
        i += 1
    clip = soft_clip + hard_clip
    
    return matches, mismatches, insertions, deletions, clip

def calculate_identities(read):
    """
    Calculate three types of identity based on CIGAR string and NM tag
    """
    if not read.has_tag('NM'):
        return None, None, None
    
    nm_tag = read.get_tag('NM')
    cigar_tuples = read.cigartuples
    
    if not cigar_tuples:
        return None, None, None
    
    # Parse CIGAR operations
    matches, mismatches, insertions, deletions, clip = parse_cigar_operations(cigar_tuples)
    
    # Adjust matches based on NM tag
    # NM = mismatches + insertions + deletions
    # If we have explicit mismatch operations (X), use them
    # if mismatches > 0:
    #     matches = matches_raw
    # else:
    #     # NM includes mismatches within M operations
    #     total_indels = insertions + deletions
    #     mismatches_in_m = nm_tag - total_indels
    #     matches = matches_raw - mismatches_in_m
    #     mismatches = mismatches_in_m
    
    # Calculate read length (aligned portion)
    read_aligned_length = matches + mismatches + insertions
    
    # Calculate reference length (aligned portion)
    ref_aligned_length = matches + mismatches + deletions
    
    # 1. Query identity (fraction of read bases that match)
    query_identity = matches / read_aligned_length if read_aligned_length > 0 else 0
    
    # 2. Reference identity (fraction of reference bases covered that match)
    reference_identity = matches / ref_aligned_length if ref_aligned_length > 0 else 0
    
    # 3. Alignment identity (fraction of aligned columns that are matches)
    total_aligned_columns = matches + mismatches + insertions + deletions
    alignment_identity = matches / total_aligned_columns if total_aligned_columns > 0 else 0

    no_indel_identity = matches / (matches + mismatches) if (matches + mismatches) > 0 else 0

    return query_identity, reference_identity, alignment_identity, matches, mismatches, insertions, deletions, no_indel_identity, clip

def handle_each_contig(bam, fig, prefix, NM_file, iden_fig, q=20):
    samfile = pysam.AlignmentFile(bam, "rb")
    f = open(NM_file, 'w')
    
    # Write header
    print("read_name\tNM\tmatches\tmismatches\tinsertions\tdeletions\tquery_identity\treference_identity\talignment_identity\tno_indel_identity", file=f)
    
    NM_list = []
    query_iden_list = []
    ref_iden_list = []
    align_iden_list = []
    no_iden_list = []
    good_ide_clip = []
    
    for read in samfile:
        if read.is_unmapped:
            continue
        if read.mapping_quality < q:
            continue
        if not read.has_tag('NM'):
            continue
            
        nm_tag = read.get_tag("NM")
        NM_list.append(nm_tag)
        
        # Calculate identities
        query_id, ref_id, align_id, matches, mismatches, insertions, deletions, no_indel_identity, clip = calculate_identities(read)

        query_iden_list.append(query_id)
        ref_iden_list.append(ref_id)
        align_iden_list.append(align_id)
        no_iden_list.append(no_indel_identity)
        if no_indel_identity > 0.99 and clip > 200: 
            clip_ratio = clip/ (matches + mismatches)
            good_ide_clip.append(clip_ratio)
            print (f"{read.query_name}\t{nm_tag}\t{matches}\t{mismatches}\t{insertions}\t{deletions}\t{query_id:.4f}\t{ref_id:.4f}\t{align_id:.4f}\t{no_indel_identity:.4f}\t{clip}")

        print(f"{read.query_name}\t{nm_tag}\t{matches}\t{mismatches}\t{insertions}\t{deletions}\t{query_id:.4f}\t{ref_id:.4f}\t{align_id:.4f}\t{no_indel_identity:.4f}", file=f)

        # if len(NM_list) % 10000 == 0 and len(NM_list) > 0:
        #     print(f"Processed {len(NM_list)} reads...")
    
    f.close()
    samfile.close()
    ## plot the distribution of good_ide_clip
    plot(good_ide_clip, "../tmp/results2/clip.pdf", "clip", "Good Indel Clip Ratio")


    # Plot all three identity types
    plot_identities(query_iden_list, ref_iden_list, align_iden_list, no_iden_list, iden_fig, prefix)
    # plot(NM_list, fig, prefix, "NM values")

def plot_identities(query_iden, ref_iden, align_iden, no_iden, fig, prefix):
    """Plot all four types of identity"""
    fig_obj, axes = plt.subplots(2, 2, figsize=(18, 10))
    axes = axes.flatten()
    identities = [
        (query_iden, "Query Identity", "Fraction of read bases that match"),
        (ref_iden, "Reference Identity", "Fraction of reference bases that match"), 
        (align_iden, "Alignment Identity", "Fraction of aligned columns that match"),
        (no_iden, "No Indel Identity", "Fraction of aligned columns that are matches (no indels)")
    ]
    
    for i, (data, title, xlabel) in enumerate(identities):
        if data:
            mean_val = round(np.mean(data), 4)
            median_val = round(np.median(data), 4)
            q1_val = round(np.percentile(data, 25), 4)
            q3_val = round(np.percentile(data, 75), 4)
            
            axes[i].hist(data, bins=50, alpha=0.7, edgecolor='black')
            axes[i].set_xlabel(xlabel)
            axes[i].set_ylabel('Read count')
            axes[i].set_title(f'{title} {prefix} \n Mean: {mean_val}, Median: {median_val}, Q1: {q1_val}, Q3: {q3_val}')
            # axes[i].set_yscale('log')
            
            print(f"{prefix} {title} - Mean: {mean_val}, Median: {median_val}, Q1: {q1_val}, Q3: {q3_val}")
        else:
            axes[i].text(0.5, 0.5, 'No data', ha='center', va='center', transform=axes[i].transAxes)
            axes[i].set_title(title)
    
    plt.tight_layout()
    plt.savefig(fig, dpi=300, bbox_inches='tight')
    plt.close()

def plot(data_list, fig, prefix, ylabel):
    """Generic plotting function"""
    if data_list:
        mean_val = round(np.mean(data_list), 2)
        median_val = round(np.median(data_list), 2)
        q1_val = round(np.percentile(data_list, 25), 2)
        q3_val = round(np.percentile(data_list, 75), 2)
        print(f"{prefix} {ylabel} - Mean: {mean_val}, Median: {median_val}, Q1: {q1_val}, Q3: {q3_val}")
    else:
        print(f"No data found for {ylabel}")
        mean_val, median_val, q1_val, q3_val = 0, 0, 0, 0

    plt.figure(figsize=(8, 6))
    sns.histplot(data_list, bins=50, kde=False)
    plt.xlabel(ylabel)
    plt.ylabel('Read count')
    plt.title(f'{prefix}\nMean: {mean_val}, Median: {median_val}, Q1: {q1_val}, Q3: {q3_val}')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(fig, dpi=300, bbox_inches='tight')
    plt.close()

# Main execution
# prefix = "ERR12723529_mice"
prefix = "infant_1"

# all_dir = "/home/shuaiw/borg/paper/run2/"
# for my_dir in os.listdir(all_dir):
#     prefix = my_dir
#     print (prefix)
bam = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.bam"
fig = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.NM.pdf"
iden_fig = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.identities.pdf"
NM_file = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.analysis.tsv"
# depth_file = f"/home/shuaiw/borg/paper/run2/{prefix}_methylation2/mean_depth.csv"

if os.path.exists(bam):
    print(f"Processing {bam}")
    handle_each_contig(bam, fig, prefix, NM_file, iden_fig)
else:
    print(f"BAM file not found: {bam}")

# prefix = "ERR12723529_mice"


# # 

# # for my_dir in os.listdir(all_dir):
#     # prefix = my_dir
#     # print (prefix)
# bam = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.bam"
# fig = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.NM.pdf"
# iden_fig = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.iden.pdf"
# NM_file = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.align.NM.csv"

# ## check if bam exists
# if os.path.exists(bam):
#     print(f"Processing {bam}")
#     handle_each_contig(bam, fig, prefix,NM_file,iden_fig)
#     # read_NM(fig, prefix, NM_file, iden_fig)