import os
import sys
import pandas as pd

def summary_host(host_dir_list, total_summary_file):
    merged_summary = pd.DataFrame()
    for host_dir in host_dir_list:
        print (f"processing {host_dir}")
        ## check if the directory exists
        if not os.path.exists(host_dir):
            print (f"directory {host_dir} does not exist")
            sys.exit(1)
        data = []
        for file in os.listdir(host_dir):
            if file.endswith(".host_prediction.csv"):
                plasmid_name = file.split(".")[0]
                host_prediction = os.path.join(host_dir, file)
                df = pd.read_csv(host_prediction)
                ## extract the first row
                ## add new column for plasmid_name at the start
                if len(df) > 0:
                    df.insert(0, 'plasmid', plasmid_name)
                    data.append(df.iloc[0])
        data = pd.DataFrame(data)
        ## sort by final_score
        ## merge the data to merged_summary, if the plasmid is already in merged_summary, then take the one with the higher final_score
        data = data.sort_values(by='final_score', ascending=False)
        merged_summary = pd.concat([merged_summary, data], ignore_index=True)
        print (f"merged_summary shape: {merged_summary.shape}")
        merged_summary = merged_summary.sort_values(by='final_score', ascending=False)
        ## drop the duplicates
        merged_summary = merged_summary.drop_duplicates(subset=['plasmid'], keep='first')
        print (f"merged_summary shape after drop duplicates: {merged_summary.shape}")
        ## sort by final_score
    merged_summary = merged_summary.sort_values(by='final_score', ascending=False)
    merged_summary.to_csv(total_summary_file, index = False)


if __name__ == "__main__":
    # host_dir_list = ["/home/shuaiw/methylation/data/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/hosts", \
    #                  "/home/shuaiw/methylation/data/borg/pengfan/RuReacBro_20230708_26_72h_NC_r4_LR_bin/hosts"]
    host_dir_list = ["/home/shuaiw/methylation/data/borg/pengfan/RuReacBro_20230708_11_72h_20_bin/hosts"]
    total_summary_file = "/home/shuaiw/methylation/data/borg/pengfan/total_summary.csv"
    summary_host(host_dir_list, total_summary_file)
    # main()