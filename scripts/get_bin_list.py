## given a bin folder, output each contig in each bin


import os
import pandas as pd
from Bio import SeqIO
import sys



def read_bins(bin_dir, bin_list_file):
    h = open(bin_list_file, 'w')
    for each_bin in os.listdir(bin_dir):
        bin_path = os.path.join(bin_dir, each_bin)
        if os.path.isfile(bin_path):
            bin_name = each_bin.split(".")[0]
            with open(bin_path, 'r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    print (f"{str(record.id)}\t{bin_name}", file=h)
        else:
            print(f"Skipping {bin_path}, not a file.")
    h.close()
if __name__ == "__main__":
    # Example usage:
    # bin_dir = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_new/bins"
    # bin_list_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_new/bin_list.txt"
    
    if len(sys.argv) != 3:
        print("Usage: python get_bin_list.py <bin_dir> <bin_list_file>")
        sys.exit(1)
    bin_dir = sys.argv[1]
    bin_list_file = sys.argv[2]
    read_bins(bin_dir, bin_list_file)