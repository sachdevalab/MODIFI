import pandas as pd



def get_edge():
    """
    Get the edge data from the host summary file.
    """
    # Read the host summary file
    host_sum = pd.read_csv(host_sum_file)
    
    # Extract the edge data
    edge = host_sum[['MGE', 'host', 'final_score']]
    
    # Save the edge data to a csv file
    edge.to_csv(edge_file, index=False)

    ## get node file
    node_data = []
    for index, row in host_sum.iterrows():
        node_data.append([row['MGE'], row['MGE'], "MGE"])
        node_data.append([row['host'], row['host'], "host"])
    node = pd.DataFrame(node_data, columns=['id', 'label', 'type'])
    node.to_csv(node_file, index=False)


if __name__ == "__main__":  
    host_sum_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/host_summary.csv"
    edge_file = "/home/shuaiw/borg/paper/network/RuReacBro_20230708_11_72h_20_bin_edge.csv"
    node_file = "/home/shuaiw/borg/paper/network/RuReacBro_20230708_11_72h_20_bin_node.csv"
    get_edge()  