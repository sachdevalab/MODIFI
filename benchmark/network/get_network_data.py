import pandas as pd



def get_edge():
    """
    Get the edge data from the host summary file.
    """
    # Read the host summary file
    host_sum = pd.read_csv(host_sum_file)
    

    ## get node file
    node_data = []
    edge_data = []
    for index, row in host_sum.iterrows():
        node_data.append([row['MGE'], row['MGE'], "MGE"])
        node_data.append([row['host'], row['host'], "host"])
        edge_data.append([row['MGE'], row['host'], row['final_score'], "Undirected"])
    node = pd.DataFrame(node_data, columns=['Id', 'Label', 'Type'])
    edge = pd.DataFrame(edge_data, columns=['Source', 'Target', 'Weight', 'Type'])
    node.to_csv(node_file, index=False)
    edge.to_csv(edge_file, index=False)


if __name__ == "__main__":  
    host_sum_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/host_summary.csv"
    edge_file = "/home/shuaiw/borg/paper/network/RuReacBro_20230708_11_72h_20_bin_edge.csv"
    node_file = "/home/shuaiw/borg/paper/network/RuReacBro_20230708_11_72h_20_bin_node.csv"
    get_edge()  