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

from repeat_count import  read_ref


def count_motif_composition(motif_seq):
    """Count GC and AT content in a motif sequence
    
    Args:
        motif_seq: DNA sequence string (only considers A, T, C, G)
        
    Returns:
        dict: Dictionary with counts and percentages
    """
    motif_seq = motif_seq.upper()
    
    # Count bases (only A, T, C, G)
    g_count = motif_seq.count('G')
    c_count = motif_seq.count('C')
    a_count = motif_seq.count('A')
    t_count = motif_seq.count('T')
    
    gc_count = g_count + c_count
    at_count = a_count + t_count
    
    # Total length based only on A, T, C, G
    length = gc_count + at_count
    
    # Calculate percentages
    gc_percent = (gc_count / length * 100) if length > 0 else 0
    at_percent = (at_count / length * 100) if length > 0 else 0
    
    return {
        'motif': motif_seq,
        'length': length,
        'G': g_count,
        'C': c_count,
        'A': a_count,
        'T': t_count,
        'GC_count': gc_count,
        'AT_count': at_count,
        'GC_percent': gc_percent,
        'AT_percent': at_percent
    }



def analyze_motif_gc_relationship(REF, host_motifs_file, window_size=5000, output_dir="../../tmp/figures/borg_fig"):
    ## please count the host-targeted methylation site density vs. GC content. 
    # Im thinking that methylation may be what motivates Borgs (etc) to have low (and for betaZs, very low) GC content..


    """Analyze the relationship between host motif occurrence density and GC content"""
    from Bio.SeqUtils import gc_fraction
    import csv
    
    # Read host motifs
    host_motifs = []
    with open(host_motifs_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            motif_str = row['motifString']
            center_pos = int(row['centerPos'])
            host_motifs.append({
                'motif': motif_str,
                'centerPos': center_pos,
                'fraction': float(row['fraction']),
                'nDetected': int(row['nDetected']),
                'nGenome': int(row['nGenome'])
            })
    
    print(f"\n{'='*60}")
    print(f"MOTIF DENSITY vs GC CONTENT ANALYSIS")
    print(f"{'='*60}")
    print(f"Analyzing {len(host_motifs)} host motifs with window size {window_size} bp\n")
    
    # Display motif composition
    print("Host Motif Composition:")
    print(f"{'Motif':<20} {'Length':<8} {'GC%':<8} {'AT%':<8} {'G':<4} {'C':<4} {'A':<4} {'T':<4}")
    print("-" * 70)
    for motif_info in host_motifs:
        comp = count_motif_composition(motif_info['motif'])
        print(f"{comp['motif']:<20} {comp['length']:<8} {comp['GC_percent']:<8.1f} "
              f"{comp['AT_percent']:<8.1f} {comp['G']:<4} {comp['C']:<4} "
              f"{comp['A']:<4} {comp['T']:<4}")
    print()
    
    # Create subplots
    n_motifs = len(host_motifs)
    n_cols = min(3, n_motifs)
    n_rows = (n_motifs + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
    if n_motifs == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_rows > 1 else axes
    
    # Analyze each motif
    for idx, motif_info in enumerate(host_motifs):
        motif = motif_info['motif']
        center_pos = motif_info['centerPos']
        
        print(f"Analyzing motif: {motif}")
        
        # Get motif positions for each contig
        all_positions = []
        all_gc = []
        
        for contig_id, seq in REF.items():
            seq_len = len(seq)
            if seq_len < window_size:
                continue
            
            # Find all motif occurrences
            motif_positions = []
            for site in nt_search(str(seq), motif)[1:]:
                motif_positions.append(site + center_pos)
            
            # Also search reverse complement
            for site in nt_search(str(seq), Seq(motif).reverse_complement())[1:]:
                rev_center = len(motif) - center_pos + 1
                motif_positions.append(site + rev_center)
            
            # Calculate motif density and GC content in sliding windows
            for start in range(0, seq_len - window_size, window_size // 2):
                end = start + window_size
                window_seq = seq[start:end]
                
                # Count motifs in this window
                motif_count = sum(1 for pos in motif_positions if start <= pos < end)
                motif_density = (motif_count / window_size) * 1000  # per kb
                
                # Calculate GC content
                gc_content = gc_fraction(window_seq) * 100
                
                all_positions.append(motif_density)
                all_gc.append(gc_content)
        
        if len(all_positions) == 0:
            print(f"  No data for motif {motif}")
            continue
        
        # Calculate correlation
        if len(all_positions) > 1:
            correlation, pvalue = pearsonr(all_positions, all_gc)
            print(f"  Pearson correlation: {correlation:.4f}, p-value: {pvalue:.4e}")
            print(f"  Mean motif density: {np.mean(all_positions):.4f} per kb")
            print(f"  Mean GC content: {np.mean(all_gc):.2f}%\n")
        
        # Plot
        ax = axes[idx]
        ax2 = ax.twinx()
        
        # Sort by GC content for better visualization
        sorted_indices = np.argsort(all_gc)
        sorted_gc = np.array(all_gc)[sorted_indices]
        sorted_density = np.array(all_positions)[sorted_indices]
        
        # Smooth the data for plotting
        window = min(50, len(sorted_gc) // 10)
        if window > 0:
            smooth_gc = np.convolve(sorted_gc, np.ones(window)/window, mode='valid')
            smooth_density = np.convolve(sorted_density, np.ones(window)/window, mode='valid')
        else:
            smooth_gc = sorted_gc
            smooth_density = sorted_density
        
        line1 = ax.plot(smooth_gc, smooth_density, 'b-', linewidth=2, label='Motif density')
        ax2.plot(smooth_gc, smooth_gc, 'r-', linewidth=2, label='GC content')
        
        ax.set_xlabel('GC content (%)', fontsize=10)
        ax.set_ylabel('Motif density (per kb)', color='b', fontsize=10)
        ax.tick_params(axis='y', labelcolor='b')
        ax2.set_ylabel('GC content (%)', color='r', fontsize=10)
        ax2.tick_params(axis='y', labelcolor='r')
        
        if len(all_positions) > 1:
            title = f"{motif}\nr={correlation:.3f}, p={pvalue:.2e}"
        else:
            title = motif
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    # Hide unused subplots
    for idx in range(n_motifs, len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'motif_density_vs_gc.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.close()
    
    print(f"{'='*60}\n")



if __name__ == "__main__":
    borg_ref = "/home/shuaiw/borg/paper/borg_data/batch_export2/BLACK_Borg-presumed-host-methylation_sites_BLACK-SR-VP_26_10_2019_C_40cm_scaffold_23_FINAL_IR.fasta"
    # borg_ref = "/home/shuaiw/borg/paper/borg_data/batch_export2/black_borgs/SR-VP_07_25_2022_A1_115cm_PACBIO-HIFI_Black_Borg_32_04.contigs.fa"
    host_motifs_file = "/home/shuaiw/borg/paper/gg_run3/soil_60/soil_60_methylation4/motifs/SR-VP_07_25_2022_A1_60cm_PACBIO-HIFI_METAMDBG_551173_L.motifs.csv"
    output_dir="../../tmp/figures/borg_fig"
    REF = read_ref(borg_ref)


    if os.path.exists(host_motifs_file):
        analyze_motif_gc_relationship(REF, host_motifs_file, window_size=10000, 
                                    output_dir=output_dir)
    else:
        print(f"\nHost motifs file not found: {host_motifs_file}")
        print("Skipping motif-GC relationship analysis")