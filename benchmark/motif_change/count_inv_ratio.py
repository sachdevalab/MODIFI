import os
import pandas as pd
import pysam


def parse_region(region_str):
    """Parse region string like 'contig:start-end' into components
    
    Args:
        region_str: String in format 'contig:start-end'
        
    Returns:
        tuple: (contig, start, end)
    """
    contig, coords = region_str.split(':')
    start, end = map(int, coords.split('-'))
    return contig, start, end


def index_bam_if_needed(bam_file):
    """Create BAM index if it doesn't exist
    
    Args:
        bam_file: Path to BAM file
    """
    if not os.path.exists(f"{bam_file}.bai"):
        print(f"  Creating index for {os.path.basename(bam_file)}...")
        os.system(f"samtools index {bam_file}")


def count_boundary_crossing_reads(bam_file, contig, start, end):
    """Count reads that cross region boundaries with quality filters
    
    Args:
        bam_file: Path to BAM file
        contig: Contig/chromosome name
        start: Region start position
        end: Region end position
        boundary_window: Window size (±bp) around boundaries to consider crossing
        max_nm: Maximum edit distance (NM tag) allowed
        
    Returns:
        int: Number of reads crossing boundaries
    """
    if not os.path.exists(bam_file):
        print(f"  Warning: {bam_file} not found")
        return 0
    
    crossing_reads = 0
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Fetch reads in expanded region to catch boundary-crossing reads
        for read in bam.fetch(contig, max(0, start - 100), end + 100):
            # Skip unmapped, secondary, and supplementary alignments
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            

            

            
            if read.reference_start < start - 20 and read.reference_end > end + 20:
                crossing_reads += 1
    
    return crossing_reads


def analyze_inversion_support(raw_region, inverted_region, raw_align_dir, inv_align_dir, 
                               weeks, samples):
    """Analyze read support for raw vs inverted regions across time points
    
    Args:
        raw_region: Region string for raw alignment
        inverted_region: Region string for inverted alignment
        raw_align_dir: Directory containing raw alignments
        inv_align_dir: Directory containing inverted alignments
        weeks: List of time point labels
        samples: List of sample names corresponding to weeks
        boundary_window: Window size for boundary crossing detection
        max_nm: Maximum edit distance filter
        
    Returns:
        list: Results dictionary for each time point
    """
    # Parse regions
    raw_contig, raw_start, raw_end = parse_region(raw_region)
    inv_contig, inv_start, inv_end = parse_region(inverted_region)
    
    print(f"{'='*70}")
    print(f"INVERSION SUPPORT ANALYSIS")
    print(f"{'='*70}")
    print(f"Raw region:      {raw_contig}:{raw_start}-{raw_end}")
    print(f"Inverted region: {inv_contig}:{inv_start}-{inv_end}")

    print(f"{'='*70}\n")
    
    results = []
    
    for week, sample in zip(weeks, samples):
        print(f"Processing {week} ({sample})...")
        
        raw_file = f"{raw_align_dir}/{week}.region.bam"
        inverted_file = f"{inv_align_dir}/{week}.region.bam"
        
        # Create indices if needed
        index_bam_if_needed(raw_file)
        index_bam_if_needed(inverted_file)
        
        # Count crossing reads
        raw_crossing = count_boundary_crossing_reads(
            raw_file, raw_contig, raw_start, raw_end
        )
        
        inv_crossing = count_boundary_crossing_reads(
            inverted_file, inv_contig, inv_start, inv_end
        )
        
        # Calculate ratio
        ratio = inv_crossing / (inv_crossing + raw_crossing )
        
        print(f"  Raw crossing reads:      {raw_crossing:>6}")
        print(f"  Inverted crossing reads: {inv_crossing:>6}")
        print(f"  Ratio (inv/raw):         {ratio:>6.4f}")
        print()
        
        results.append({
            'sample': sample,
            'week': week,
            'raw_crossing': raw_crossing,
            'inv_crossing': inv_crossing,
            'ratio': ratio
        })
    
    return results


def print_summary(results):
    """Print summary table of results
    
    Args:
        results: List of result dictionaries
    """
    print(f"\n{'='*80}")
    print(f"SUMMARY")
    print(f"{'='*80}")
    print(f"{'Sample':<12} {'Week':<12} {'Raw':<15} {'Inverted':<15} {'Ratio':<15}")
    print(f"{'-'*80}")
    
    for r in results:
        print(f"{r['sample']:<12} {r['week']:<12} {r['raw_crossing']:<15} {r['inv_crossing']:<15} {r['ratio']:<15.4f}")
    
    print(f"{'='*80}\n")


def main():
    """Main execution function"""
    # Configuration
    raw_region = "infant_2_3_C:2894667-2897926"
    inverted_region = "infant_14_31_C:2323856-2327115"
    raw_align = "/home/shuaiw/borg/paper/motif_change/new_alignment/"
    inverted_align = "/home/shuaiw/borg/paper/motif_change/new_alignment2/"
    weeks = ["week2", "week3", "week4", "week5"]
    samples = ["DOL14","DOL21", "DOL28", "DOL35"]
    
    # Run analysis
    results = analyze_inversion_support(
        raw_region=raw_region,
        inverted_region=inverted_region,
        raw_align_dir=raw_align,
        inv_align_dir=inverted_align,
        weeks=weeks,
        samples=samples
    )
    
    # Print summary
    print_summary(results)
    ## save the results as df
    df = pd.DataFrame(results)
    df.to_csv("../../tmp/figures/inversion/inversion_ratio_results.csv", index=False)


if __name__ == "__main__":
    main()
