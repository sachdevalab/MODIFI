def read_depth(depth_file):
    """
    Read the depth file and return a dictionary mapping borg names to their depths.
    """
    borg_depth = {}
    with open(depth_file, 'r') as f:
        f.readline()
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 2:
                borg_name, depth = parts
                borg_depth[borg_name] = float(depth)
    return borg_depth

def count_borg_dp(depth_file, borg_list):
    borg_depth = read_depth(depth_file)
    count = 0
    with open(borg_list, 'r') as f:
        f.readline()
        for line in f:
            borg_name = line.strip().split()[0]
            if borg_name in borg_depth:
                print (f"{borg_name}: {borg_depth[borg_name]}")


depth_file = "/home/shuaiw/borg/bench/soil_borg/mean_depth.csv"
# depth_file = "/home/shuaiw/borg/bench/soil/run1/mean_depth.csv"
# depth_file = "/home/shuaiw/borg/bench/soil/m84039_230626_214113_s4.hifi_reads.bc2022/mean_depth.csv"
# borg_list = "/home/shuaiw/borg/bench/soil_borg/borg.list"
borg_list = "/home/shuaiw/borg/bench/soil_borg/host.list"
count_borg_dp(depth_file, borg_list)