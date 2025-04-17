"""  
find {input.ref_dir} -name "*.fa" > {work_dir}/contigs_list.txt"
"""

import os
import sys




def get_contig_list(ipd_dir, ref_dir, list_file):
    f = open(list_file, "w")
    for file in os.listdir(ipd_dir):
        if file.endswith(".ipd1.csv"):
            ## sample name
            sample = file.split(".")[0]
            contig = os.path.join(ref_dir, sample + ".fa")
            if not os.path.exists(contig):
                print (f"Contig {contig} not found")
                continue
            else:
                f.write(contig + "\n")
    f.close()
    
if __name__ == "__main__":
    work_dir = sys.argv[1]

    ref_dir = os.path.join(work_dir, "contigs")
    ipd_dir = os.path.join(work_dir, "ipd")
    list_file = os.path.join(work_dir, "contigs_list.txt")

    get_contig_list(ipd_dir, ref_dir, list_file)



