import os



## remove repeat contigs and only keep these with length > 10000
def filter_all_borgs(raw_merge, filter_merge):
    cmd = f"seqkit rmdup -s -i {raw_merge} | seqkit seq -m 10000 > {filter_merge}"
    os.system(cmd)

if __name__ == "__main__":
    raw_merge = "/home/shuaiw/borg/paper/borg_data/all_borgs.raw.fa"
    filter_merge = "/home/shuaiw/borg/paper/borg_data/all_borgs.fa"
    filter_all_borgs(raw_merge, filter_merge)