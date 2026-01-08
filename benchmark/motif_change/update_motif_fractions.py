#!/usr/bin/env python3

import pandas as pd

def update_motif_fractions(metadata_file, fraction_file, output_file):
    """
    Update motif fraction values in metadata file using values from fraction file.
    
    Args:
        metadata_file: Path to 180_4.csv (has sample, type, anno, contig, and old motif fractions)
        fraction_file: Path to 161_4.csv (has contig, motifString, fraction in long format)
        output_file: Path to save updated CSV
    """
    # Read fraction data into dictionary: {(contig, motifString): fraction}
    print(f"Reading fraction data from {fraction_file}...")
    fraction_df = pd.read_csv(fraction_file)
    print(f"Fraction data shape: {fraction_df.shape}")
    
    fraction_dict = {}
    for _, row in fraction_df.iterrows():
        contig = row['contig']
        motif = row['motifString']
        fraction = row['fraction']
        fraction_dict[(contig, motif)] = fraction
    
    print(f"Loaded {len(fraction_dict)} contig-motif pairs into dictionary")
    
    # Read metadata file
    print(f"\nReading metadata from {metadata_file}...")
    metadata_df = pd.read_csv(metadata_file)
    print(f"Metadata shape: {metadata_df.shape}")
    print(f"Metadata columns: {metadata_df.columns.tolist()}")
    
    # Get metadata columns
    meta_columns = ['sample', 'type', 'anno', 'contig']
    motif_columns = [col for col in metadata_df.columns if col not in meta_columns]
    print(f"Found {len(motif_columns)} motif columns")
    
    # Update fraction values
    print("\nUpdating motif fractions...")
    updated_count = 0
    not_found_count = 0
    
    for idx, row in metadata_df.iterrows():
        contig = row['contig']
        
        # Fix anno: change P2 to P1 where type is plasmid_1
        if 'plasmid_1' in row['type'] and row['anno'].endswith('_P2'):
            metadata_df.at[idx, 'anno'] = row['anno'].replace('_P2', '_P1')
        
        # Update each motif column
        for motif in motif_columns:
            key = (contig, motif)
            if key in fraction_dict:
                metadata_df.at[idx, motif] = fraction_dict[key]
                updated_count += 1
            else:
                # Set to 0 if motif not found for this contig
                metadata_df.at[idx, motif] = 0.0
                not_found_count += 1
    
    print(f"Updated {updated_count} values")
    print(f"Set {not_found_count} missing values to 0")
    
    # Save to output file
    metadata_df.to_csv(output_file, index=False)
    print(f"\nUpdated data saved to: {output_file}")
    
    # Show sample of updated data
    print("\nSample of updated data:")
    print(metadata_df.head())
    
    # Show comparison for first contig
    if len(metadata_df) > 0:
        first_contig = metadata_df['contig'].iloc[0]
        print(f"\nExample values for contig: {first_contig}")
        first_row = metadata_df[metadata_df['contig'] == first_contig].iloc[0]
        for col in motif_columns[:5]:  # Show first 5 motifs
            print(f"  {col}: {first_row[col]:.6f}")

if __name__ == "__main__":
    metadata_file = "180_4.csv"
    fraction_file = "161_4.csv"
    output_file = "180_4_updated.csv"
    
    update_motif_fractions(metadata_file, fraction_file, output_file)
    print("\nDone! You can now use 180_4_updated.csv with the R script.")
