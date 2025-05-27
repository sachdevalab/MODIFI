import pandas as pd
import networkx as nx



def get_edge(cutoff=0.45):
    """
    Get the edge data from the host summary file.
    """
    # Read the host summary file
    host_sum = pd.read_csv(host_sum_file)
    
    G = nx.Graph()
    ## get node file
    for index, row in host_sum.iterrows():
        if row['final_score'] < cutoff:
            continue

        G.add_node(row['MGE'], label=row['MGE'], type="MGE")
        G.add_node(row['host'], label=row['host'], type="host")
        G.add_edge(row['MGE'], row['host'], weight=row['final_score'])

    nx.write_gexf(G,gexf)


if __name__ == "__main__":  
    # host_sum_file = "/home/shuaiw/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/host_summary.csv"
    # gexf = "/home/shuaiw/borg/paper/network/RuReacBro_20230708_11_72h_20_bin.gexf"

    host_sum_file = "/home/shuaiw/borg/bench/soil/run1/host_summary.csv"
    gexf = "/home/shuaiw/borg/paper/network/soil_run1_bin.gexf"
    get_edge()  