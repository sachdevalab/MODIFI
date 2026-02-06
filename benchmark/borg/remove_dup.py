import os
from Bio import SeqIO

folder = "/home/shuaiw/borg/paper/curated_genome/"
curated_genome_list = f"{folder}/curated_genome.list"
output_folder = f"{folder}/unique"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Read the curated genome list and collect all duplicate genome IDs
duplicate_ids = set()
with open(curated_genome_list, 'r') as f:
    for line in f:
        parts = line.strip().split()

        if len(parts) >= 6:
            # Last column contains comma-separated genome IDs to remove
            dup_genomes = parts[-1].split(',')
            for genome_id in dup_genomes:
                duplicate_ids.add(genome_id.strip())

print(f"Found {len(duplicate_ids)} duplicate genome IDs to remove")

# Process each .fa file in the folder
fa_files = [f for f in os.listdir(folder) if f.endswith('.fa')]
print(f"Found {len(fa_files)} .fa files to process")

for fa_file in fa_files:
    input_path = os.path.join(folder, fa_file)
    output_path = os.path.join(output_folder, fa_file)
    
    # Filter sequences
    kept_records = []
    removed_count = 0
    
    for record in SeqIO.parse(input_path, "fasta"):
        # Check if the sequence ID is in the duplicate list
        if record.id not in duplicate_ids:
            kept_records.append(record)
        else:
            removed_count += 1
    
    # Write filtered sequences to output
    SeqIO.write(kept_records, output_path, "fasta")
    
    print(f"{fa_file}: kept {len(kept_records)} sequences, removed {removed_count} duplicates")

print(f"\nFiltered genomes saved to {output_folder}")