#!/usr/bin/env python3
"""
Given a BAM file, count how many reads are mapped in reverse orientation
"""

import pysam
import sys
from collections import defaultdict

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

def print_statistics(stats):
    """Print detailed statistics about read orientations"""
    
    print("\n" + "="*50)
    print("READ ORIENTATION STATISTICS")
    print("="*50)
    
    print(f"Total reads: {stats['total_reads']:,}")
    print(f"Mapped reads: {stats['mapped_reads']:,}")
    print(f"Unmapped reads: {stats['unmapped_reads']:,}")
    
    if stats['mapped_reads'] > 0:
        forward_pct = (stats['forward_mapped'] / stats['mapped_reads']) * 100
        reverse_pct = (stats['reverse_mapped'] / stats['mapped_reads']) * 100
        
        print(f"\nMapped Read Orientations:")
        print(f"  Forward strand (+): {stats['forward_mapped']:,} ({forward_pct:.1f}%)")
        print(f"  Reverse strand (-): {stats['reverse_mapped']:,} ({reverse_pct:.1f}%)")
        
        # Print per-reference statistics
        if len(stats['by_reference']) > 1:
            print(f"\nPer-reference statistics:")
            print(f"{'Reference':<30} {'Forward':<10} {'Reverse':<10} {'Total':<10} {'% Reverse':<10}")
            print("-" * 80)
            
            for ref_name in sorted(stats['by_reference'].keys()):
                ref_stats = stats['by_reference'][ref_name]
                total_ref = ref_stats['forward'] + ref_stats['reverse']
                reverse_pct_ref = (ref_stats['reverse'] / total_ref * 100) if total_ref > 0 else 0
                
                print(f"{ref_name:<30} {ref_stats['forward']:<10} {ref_stats['reverse']:<10} "
                      f"{total_ref:<10} {reverse_pct_ref:<10.1f}%")

def main():
    if len(sys.argv) != 2:
        print("Usage: python count_rev_map.py <bam_file>")
        print("Example: python count_rev_map.py reads.bam")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    
    try:
        stats = count_reverse_mappings(bam_file)
        print_statistics(stats)
        
    except FileNotFoundError:
        print(f"Error: BAM file '{bam_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing BAM file: {e}")
        sys.exit(1)

# Default BAM file for testing
sample="soil_2"
# default_bam = "/home/shuaiw/borg/paper/borg_data/borg_for/soil_1/soil_1_methylation3/bams/SRVP18_trench_6_60cm_scaf_214_117_86_FINAL.bam"
default_bam = f"/home/shuaiw/borg/paper/borg_data/borg_for/{sample}/{sample}_methylation3/bams/Methanoperedens_44_19-type_SR-VP_26_10_2019_1_100cm_part_1.bam"

if __name__ == "__main__":
    if len(sys.argv) == 1:
        # Use default BAM file if no argument provided
        print("No BAM file provided, using default...")
        sys.argv.append(default_bam)
    main()


