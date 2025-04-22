import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import List, Tuple
import sys
import os
import numpy as np

# -------------------------------
# Parameters
# -------------------------------
DEPTH_THRESHOLD = 5
MAX_GAP = 10


# -------------------------------
# Function: Segment genome by depth
# -------------------------------
def segment_by_depth(depth_file: str, depth_threshold: int, max_gap: int) -> List[Tuple[str, int, int]]:
    depth_df = pd.read_csv(depth_file, sep="\t", header=None, names=["seqname", "pos", "depth"])
    ## calculate the mean depth of the contig 
    mean_depth = round(np.mean(depth_df["depth"]), 2)
    print (f"mean depth: {mean_depth}", depth_file)
    segments = []
    current_seq = None
    start = None
    last_pos = None
    gap_count = 0

    for _, row in depth_df.iterrows():
        seq, pos, depth = row["seqname"], row["pos"], row["depth"]
        if seq != current_seq:
            if start is not None:
                segments.append((current_seq, start, last_pos))
            current_seq = seq
            start = None
            gap_count = 0

        if depth >= depth_threshold:
            if start is None:
                start = pos
            last_pos = pos
            gap_count = 0
        else:
            if start is not None:
                gap_count += 1
                if gap_count > max_gap:
                    segments.append((seq, start, last_pos))
                    start = None
                    gap_count = 0

    if start is not None:
        segments.append((current_seq, start, last_pos))

    


    return segments, mean_depth


# -------------------------------
# Function: Extract high-depth segments from FASTA
# -------------------------------
def extract_segmented_fasta(segments: List[Tuple[str, int, int]], fasta_path: str, output_path: str) -> List[Tuple[str, str, int, int]]:
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    segment_metadata = []

    with open(output_path, "w") as out_fa:
        for idx, (seq, start, end) in enumerate(segments):
            segment_name = f"{seq}_seg{idx+1}"
            segment_seq = fasta_dict[seq].seq[start-1:end]
            SeqIO.write(SeqRecord(segment_seq, id=segment_name, description=""), out_fa, "fasta")
            segment_metadata.append((segment_name, seq, start, end))

    return segment_metadata


# -------------------------------
# Function: Update GFF based on segments
# -------------------------------
def update_gff(input_gff: str, output_gff: str, segment_metadata: List[Tuple[str, str, int, int]]) -> None:
    segment_df = pd.DataFrame(segment_metadata, columns=["segment_name", "original_seq", "seg_start", "seg_end"])

    with open(input_gff) as gff_in, open(output_gff, "w") as gff_out:
        for line in gff_in:
            if line.startswith("#"):
                gff_out.write(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqname, start, end = parts[0], int(parts[3]), int(parts[4])
            matches = segment_df[
                (segment_df["original_seq"] == seqname) &
                (segment_df["seg_start"] <= start) &
                (segment_df["seg_end"] >= end)
            ]
            for _, row in matches.iterrows():
                new_seqname = row["segment_name"]
                new_start = start - row["seg_start"] + 1
                new_end = end - row["seg_start"] + 1
                parts[0] = new_seqname
                parts[3] = str(new_start)
                parts[4] = str(new_end)
                gff_out.write("\t".join(parts) + "\n")


# -------------------------------
# Main Pipeline Function
# -------------------------------
def process_depth_and_gff(depth_file: str, reference_fasta: str, input_gff: str,
                          output_fasta: str, output_gff: str,
                          depth_threshold: int = DEPTH_THRESHOLD, max_gap: int = MAX_GAP):
    segments, mean_depth = segment_by_depth(depth_file, depth_threshold, max_gap)
    segment_metadata = extract_segmented_fasta(segments, reference_fasta, output_fasta)
    update_gff(input_gff, output_gff, segment_metadata)
    return mean_depth

if __name__ == "__main__":
    # -------------------------------
    # Example usage:
    # process_depth_and_gff("depth.txt", "reference.fa", "input.gff", "segmented.fasta", "segmented.gff")

    # depth_file = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30/bams/E_coli_H10407_1.depth"
    # reference_fasta = "/home/shuaiw/methylation/data/borg/bench/zymo_new_ref_p0.05_cov1_s30/contigs/E_coli_H10407_1.fa"
    # raw_gff = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30_filter/gffs/E_coli_H10407_1.gff"
    # segment = "/home/shuaiw/methylation/data/borg/bench/test/segmented.fasta"
    # output_gff = "/home/shuaiw/methylation/data/borg/bench/test/segmented.gff"

    depth_file = sys.argv[1]
    reference_fasta = sys.argv[2]
    raw_gff = sys.argv[3]
    segment = sys.argv[4]
    output_gff = sys.argv[5]


    process_depth_and_gff(depth_file, reference_fasta, raw_gff, segment, output_gff,
                        depth_threshold=DEPTH_THRESHOLD, max_gap=MAX_GAP)