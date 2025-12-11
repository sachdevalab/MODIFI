#!/usr/bin/env python3
"""
BORG Detection Script

This script maps BORG reference sequences onto assembly contigs using minimap2
and identifies potential BORG sequences based on alignment quality and coverage.
"""

import subprocess
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'isolation'))
from sample_object import My_sample, Isolation_sample

def create_output_directory(work_dir):
    """Create output directory if it doesn't exist."""
    Path(work_dir).mkdir(parents=True, exist_ok=True)

def run_minimap2(assembly_fasta, work_dir, borg_ref, threads=10):
    """
    Run minimap2 to align BORG references against assembly contigs.
    
    Args:
        assembly_fasta (str): Path to assembly FASTA file
        work_dir (str): Working directory for output files
        borg_ref (str): Path to BORG reference FASTA file
        threads (int): Number of threads to use
        
    Returns:
        str: Path to the PAF output file
    """
    paf_output = os.path.join(work_dir, "minimap2.paf")
    
    # Check if input files exist
    for file_path in [assembly_fasta, borg_ref]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Input file not found: {file_path}")
    
    print(f"Running minimap2 alignment...")
    print(f"  Assembly: {assembly_fasta}")
    print(f"  BORG ref: {borg_ref}")
    print(f"  Output: {paf_output}")
    
    cmd = [
        "minimap2",
        "-t", str(threads),
        "-x", "asm5",  # Use assembly-to-assembly preset
        "--secondary=no",  # Only primary alignments
        borg_ref,
        assembly_fasta
    ]
    
    try:
        with open(paf_output, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, 
                                  text=True, check=True)
        print(f"✓ Minimap2 completed successfully")
        return paf_output
    except subprocess.CalledProcessError as e:
        print(f"✗ Minimap2 failed: {e}")
        print(f"Error output: {e.stderr}")
        sys.exit(1)

def parse_paf_file(paf_file, min_identity=0.8, min_coverage=0.5):
    """
    Parse PAF file and filter alignments based on quality criteria.
    
    Args:
        paf_file (str): Path to PAF file
        min_identity (float): Minimum alignment identity (0-1)
        min_coverage (float): Minimum query coverage (0-1)
        
    Returns:
        pd.DataFrame: Filtered alignments
    """
    if not os.path.exists(paf_file):
        raise FileNotFoundError(f"PAF file not found: {paf_file}")
    
    # PAF format columns
    columns = [
        'query_name', 'query_length', 'query_start', 'query_end',
        'strand', 'target_name', 'target_length', 'target_start', 'target_end',
        'matches', 'alignment_length', 'mapping_quality'
    ]
    
    try:
        # Read PAF file
        alignments = []
        with open(paf_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 12:  # Ensure minimum PAF fields
                        alignments.append(fields[:12])
        
        if not alignments:
            print("⚠ No alignments found in PAF file")
            return pd.DataFrame()
        
        df = pd.DataFrame(alignments, columns=columns)
        
        # Convert numeric columns
        numeric_cols = ['query_length', 'query_start', 'query_end', 
                       'target_length', 'target_start', 'target_end',
                       'matches', 'alignment_length', 'mapping_quality']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col])
        
        # Calculate metrics
        df['identity'] = df['matches'] / df['alignment_length']
        df['query_coverage'] = (df['query_end'] - df['query_start']) / df['query_length']
        df['target_coverage'] = (df['target_end'] - df['target_start']) / df['target_length']
        
        # Filter alignments
        initial_count = len(df)
        df_filtered = df[
            (df['identity'] >= min_identity) & 
            (df['query_coverage'] >= min_coverage)
        ]
        
        print(f"Alignments: {initial_count} total, {len(df_filtered)} passed filters")
        print(f"  Filters: identity ≥ {min_identity}, coverage ≥ {min_coverage}")
        
        return df_filtered
        
    except Exception as e:
        print(f"✗ Error parsing PAF file: {e}")
        return pd.DataFrame()

def identify_borg_contigs(df_alignments, output_file, prefix, all_dir):
    """
    Identify and save potential BORG contigs.
    
    Args:
        df_alignments (pd.DataFrame): Filtered alignments
        output_file (str): Path to output file
        
    Returns:
        list: List of BORG contig names
    """
    if df_alignments.empty:
        print("⚠ No BORG contigs identified")
        return []
    sample_obj = My_sample(prefix, all_dir )
    sample_obj.read_depth()
    # Group by target (assembly contig) and get best alignment per contig
    borg_contigs = []
    contig_summary = []
    
    for target_name, group in df_alignments.groupby('target_name'):
        best_alignment = group.loc[group['identity'].idxmax()]
        my_type = "BORG"
        if re.search("Methanoperedens", target_name):
            my_type = "HOST"
            # print ("find host, skip", target_name)
            # continue
        borg_contigs.append(target_name)
        contig_summary.append({
            'borg_ref': target_name,
            'type': my_type,
            'seq_name': best_alignment['query_name'],
            'identity': round(best_alignment['identity'], 3),
            'query_coverage': round(best_alignment['query_coverage'], 3),
            'target_coverage': round(best_alignment['target_coverage'], 3),
            'alignment_length': best_alignment['alignment_length'],
            'length': best_alignment['target_length'],
            'ctg_depth': sample_obj.depth_dict.get(best_alignment['query_name'], 'NA')
        })
    # Save summary
    summary_df = pd.DataFrame(contig_summary)
    summary_df = summary_df.sort_values('identity', ascending=False)
    summary_df.to_csv(output_file, index=False, sep = "\t")

    print(f"✓ Identified {len(borg_contigs)} potential BORG contigs")
    print(f"✓ Summary saved to: {output_file}")
    print (summary_df)
    
    return borg_contigs

def find_borg_func(assembly_fasta, work_dir, borg_ref, prefix,all_dir, ece_type, threads=10, 
                   min_identity=0.8, min_coverage=0.5):
    """
    Main function to find BORG sequences in assembly.
    
    Args:
        assembly_fasta (str): Path to assembly FASTA file
        work_dir (str): Working directory for output files
        borg_ref (str): Path to BORG reference FASTA file
        threads (int): Number of threads for minimap2
        min_identity (float): Minimum alignment identity
        min_coverage (float): Minimum query coverage
        
    Returns:
        list: List of BORG contig names
    """
    print(f"🔍 Starting BORG detection...")
    
    # Create output directory
    create_output_directory(work_dir)
    
    # Run minimap2
    paf_file = run_minimap2(assembly_fasta, work_dir, borg_ref, threads)
    
    # Parse alignments
    df_alignments = parse_paf_file(paf_file, min_identity, min_coverage)
    
    # Identify BORG contigs
    output_file = os.path.join(work_dir, f"{ece_type}_contigs_summary.tsv")
    borg_contigs = identify_borg_contigs(df_alignments, output_file, prefix, all_dir)
    
    # Save simple list
    borg_list_file = os.path.join(work_dir, f"{ece_type}_contigs.txt")
    with open(borg_list_file, 'w') as f:
        for contig in borg_contigs:
            f.write(f"{contig}\n")
    
    print(f"✓ {ece_type} contig list saved to: {borg_list_file}")
    print(f"🎉 {ece_type} detection completed!")
    
    return borg_contigs

def main():
    """Command line interface."""
    parser = argparse.ArgumentParser(
        description="Find BORG sequences in assembly using minimap2",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("assembly", help="Path to assembly FASTA file")
    parser.add_argument("output_dir", help="Output directory for results")
    parser.add_argument("--borg_ref", 
                       default="/home/shuaiw/borg/paper/borg_data/all_borgs_host_test.fa",
                       help="Path to BORG reference FASTA file")
    parser.add_argument("--prefix", 
                       help="prefix")
    parser.add_argument("--threads", type=int, default=30, 
                       help="Number of threads for minimap2")
    parser.add_argument("--min_identity", type=float, default=0.8,
                       help="Minimum alignment identity (0-1)")
    parser.add_argument("--min_coverage", type=float, default=0.5,
                       help="Minimum query coverage (0-1)")
    
    args = parser.parse_args()
    
    # Run BORG detection
    borg_contigs = find_borg_func(
        assembly_fasta=args.assembly,
        work_dir=args.output_dir,
        borg_ref=args.borg_ref,
        prefix=args.prefix,
        threads=args.threads,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage
    )
    

if __name__ == "__main__":
    # Default parameters for backward compatibility
    borg_ref = "/home/shuaiw/borg/paper/borg_data/borgs_mp_nanopore.contigs.fa"
    
    if len(sys.argv) == 1:
        # Run with default parameters if no arguments provided
        print("Running with default parameters...")
        # borg_ref = "/home/shuaiw/borg/paper/borg_data/jumbo_phage.fa"
        all_dir = "/home/shuaiw/borg/paper/run2/"
        for prefix in ["soil_1", "soil_2", "soil_s1_1","soil_s1_2","soil_s3_1","soil_s3_2","soil_s4_1","soil_s4_2"]:
            work_dir = f"/home/shuaiw/borg/paper/run2/{prefix}/borg/"
            assembly_fasta = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa"
            find_borg_func(assembly_fasta, work_dir, borg_ref, prefix, all_dir, ece_type="borg")
    else:
        # Use command line interface
        all_dir = "/home/shuaiw/borg/paper/run2/"
        main()
