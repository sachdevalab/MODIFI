input_gff = "/home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30_filter/gffs/E_coli_H10407_1.gff"
output_gff = "/home/shuaiw/methylation/data/borg/bench/test/filtered_under_30kbp.gff"

with open(input_gff) as fin, open(output_gff, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        start, end = int(parts[3]), int(parts[4])
        if end < 30000:
            fout.write(line)

