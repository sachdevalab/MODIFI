"""
extract the contigs with depth >= 5 and output to a new fasta file
"""

def read_depth():
    depth = {}
    with open(depth_file, "r") as f:
        f.readline()  # skip header
        for line in f:
            if line.startswith("#"):
                continue
            contig, mean_depth = line.strip().split(",")
            depth[contig] = float(mean_depth)
    return depth

def extract_high_cov_contigs(depth, assembly, output, min_depth=5):
    from Bio import SeqIO
    selected = set([contig for contig, dp in depth.items() if dp >= min_depth])
    with open(assembly, "r") as infile, open(output, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in selected:
                SeqIO.write(record, outfile, "fasta")

prefix = "soil_1"
depth_file = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_methylation_0.97/mean_depth.csv"
assembly = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}.hifiasm.p_ctg.rename.fa"
high_dp_assembly = f"/home/shuaiw/borg/paper/run2/{prefix}/{prefix}_high_dp.fa"

if __name__ == "__main__":
    depth = read_depth()
    extract_high_cov_contigs(depth, assembly, high_dp_assembly, min_depth=5)