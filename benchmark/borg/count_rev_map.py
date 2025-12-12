#!/usr/bin/env python3
"""
Given a BAM file, count how many reads are mapped in reverse orientation
"""

import pysam
import sys
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def count_reverse_mappings(bam_file):
    """
    Count forward and reverse mapped reads in a BAM file
    
    Args:
        bam_file (str): Path to BAM file
        
    Returns:
        dict: Statistics about read orientations
    """
    
    # Open BAM file
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'forward_mapped': 0,
        'reverse_mapped': 0,
        'unmapped_reads': 0,
        'by_reference': defaultdict(lambda: {'forward': 0, 'reverse': 0})
    }
    
    print(f"Processing BAM file: {bam_file}")
    
    # Count reads
    for read in samfile.fetch():
        stats['total_reads'] += 1
        
        if read.is_unmapped:
            stats['unmapped_reads'] += 1
            continue
            
        stats['mapped_reads'] += 1
        
        # Check if read is reverse complemented
        if not read.is_reverse:
            stats['reverse_mapped'] += 1
            stats['by_reference'][read.reference_name]['reverse'] += 1
        else:
            stats['forward_mapped'] += 1
            stats['by_reference'][read.reference_name]['forward'] += 1
    
    samfile.close()
    return stats

def print_statistics(stats, bam_file):
    """Print detailed statistics about read orientations"""
    
    print("\n" + "="*50)
    print("READ ORIENTATION STATISTICS")
    print("="*50)
    
    print(f"Total reads: {stats['total_reads']:,}")
    print(f"Mapped reads: {stats['mapped_reads']:,}")
    print(f"Unmapped reads: {stats['unmapped_reads']:,}")
    
    data= []
    if stats['mapped_reads'] > 1000:
        forward_pct = (stats['forward_mapped'] / stats['mapped_reads']) * 100
        reverse_pct = (stats['reverse_mapped'] / stats['mapped_reads']) * 100
        
        # print(f"\nMapped Read Orientations:")
        # print(f"  Forward strand (+): {stats['forward_mapped']:,} ({forward_pct:.1f}%)")
        # print(f"  Reverse strand (-): {stats['reverse_mapped']:,} ({reverse_pct:.1f}%)")
        
        # Print per-reference statistics
        if len(stats['by_reference']) > 0:
            # print(f"\nPer-reference statistics:")
            # print(f"{'Reference':<30} {'Forward':<10} {'Reverse':<10} {'Total':<10} {'% Reverse':<10}")
            # print("-" * 80)
            
            for ref_name in sorted(stats['by_reference'].keys()):
                ref_stats = stats['by_reference'][ref_name]
                total_ref = ref_stats['forward'] + ref_stats['reverse']
                reverse_pct_ref = (ref_stats['reverse'] / total_ref * 100) if total_ref > 0 else 0
                
                # print(f"{ref_name:<30} {ref_stats['forward']:<10} {ref_stats['reverse']:<10} "
                #       f"{total_ref:<10} {reverse_pct_ref:<10.1f}%")
                data.append([ref_name, ref_stats['forward'], ref_stats['reverse'], total_ref, reverse_pct_ref])
    
    # Plot distribution of reverse percentage across references
    if data:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
        
        df = pd.DataFrame(data, columns=['reference', 'forward', 'reverse', 'total', 'reverse_pct'])
        
        # Create distribution plot
        plt.figure(figsize=(12, 6))
        
        # Subplot 1: Histogram of reverse percentages
        plt.subplot(1, 2, 1)
        plt.hist(df['reverse_pct'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        plt.xlabel('Reverse Percentage (%)')
        plt.ylabel('Number of References')
        plt.title('Distribution of Reverse Read Percentages')
        plt.grid(True, alpha=0.3)
        
        # Add statistics on the plot
        mean_pct = df['reverse_pct'].mean()
        median_pct = df['reverse_pct'].median()
        plt.axvline(mean_pct, color='red', linestyle='--', label=f'Mean: {mean_pct:.1f}%')
        plt.axvline(median_pct, color='orange', linestyle='--', label=f'Median: {median_pct:.1f}%')
        plt.legend()
        
        # Subplot 2: Box plot
        plt.subplot(1, 2, 2)
        plt.boxplot(df['reverse_pct'])
        plt.ylabel('Reverse Percentage (%)')
        plt.title('Box Plot of Reverse Read Percentages')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save the plot
        output_file = bam_file.replace('.bam', '_reverse_distribution.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nDistribution plot saved to: {output_file}")
        
        # Print summary statistics
        print(f"\nSummary Statistics:")
        print(f"Mean reverse %: {mean_pct:.2f}%")
        print(f"Median reverse %: {median_pct:.2f}%")
        print(f"Std deviation: {df['reverse_pct'].std():.2f}%")
        print(f"Min reverse %: {df['reverse_pct'].min():.2f}%")
        print(f"Max reverse %: {df['reverse_pct'].max():.2f}%")
        
        plt.show()
    ## print top ten genomes with highest or lowest reverse_pct_ref
        top10_highest = df.nlargest(10, 'reverse_pct')
        top10_lowest = df.nsmallest(10, 'reverse_pct')
        print("\nTop 10 genomes with highest reverse %:")
        print(top10_highest[['reference', 'reverse_pct']])
        print("\nTop 10 genomes with lowest reverse %:")
        print(top10_lowest[['reference', 'reverse_pct']])


def main(bam_file):
    if not bam_file:
        print("Usage: python count_rev_map.py <bam_file>")
        print("Example: python count_rev_map.py reads.bam")
        sys.exit(1)
    
    
    try:
        stats = count_reverse_mappings(bam_file)
        print_statistics(stats, bam_file)
        
    except FileNotFoundError:
        print(f"Error: BAM file '{bam_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing BAM file: {e}")
        sys.exit(1)

# # Default BAM file for testing
# sample="soil_60"
# # default_bam = "/home/shuaiw/borg/paper/borg_data/borg_for/soil_1/soil_1_methylation3/bams/SRVP18_trench_6_60cm_scaf_214_117_86_FINAL.bam"
# default_bam = f"/home/shuaiw/borg/paper/gg_run/{sample}/{sample}_methylation3/align_bam/aligned.bam"
# # default_bam = f"/home/shuaiw/borg/paper/borg_data/borg_for/{sample}/{sample}_methylation3/align_bam/aligned.bam"

if __name__ == "__main__":
    for sample in ["soil_60", "soil_80","soil_90","soil_100","soil_110","soil_115","soil_1","soil_2"]:
        default_bam = f"/home/shuaiw/borg/paper/gg_run/{sample}/{sample}_methylation3/align_bam/aligned.bam"
        main(default_bam)


