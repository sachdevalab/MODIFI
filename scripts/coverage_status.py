import pysam
import numpy as np
import argparse

def compute_coverage_stats(bam_file, output_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    contig_stats = {}

    for contig in bam.references:
        contig_length = bam.get_reference_length(contig)
        coverage = np.zeros(contig_length, dtype=int)

        for pileupcolumn in bam.pileup(contig, truncate=True):
            pos = pileupcolumn.reference_pos
            coverage[pos] = pileupcolumn.nsegments

        mean = np.mean(coverage)
        var = np.var(coverage)

        contig_stats[contig] = (mean, var)

    # Write to output file
    with open(output_file, "w") as f:
        f.write("contig\tmean\tvariance\n")
        for contig, (mean, var) in contig_stats.items():
            f.write(f"{contig}\t{mean:.4f}\t{var:.4f}\n")

    print(f"Coverage stats written to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute mean and variance of base coverage per contig.")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output stats file")
    args = parser.parse_args()

    compute_coverage_stats(args.input, args.output)

