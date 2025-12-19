"""
for file in /home/shuaiw/borg/paper/gg_run2/soil_*/soil_*_methylation4/host_summary.csv
read them as pd dataframes, concatenate them together
output the combined dataframe to ./soil_all_hosts.summary.csv
"""


import os
import pandas as pd

def summarize_hosts(results_dir, output_file):
    df_list = []
    for folder in os.listdir(results_dir):
        prefix = folder
        summary_file = os.path.join(results_dir, folder, f"{prefix}_methylation4/host_summary.csv")
        if os.path.exists(summary_file) and os.path.getsize(summary_file) > 0:
            ## line count > 1
            line_count = sum(1 for line in open(summary_file))
            if line_count < 2:
                continue
            df = pd.read_csv(summary_file)
            df['sample'] = prefix
            df_list.append(df)
    if df_list:
        combined_df = pd.concat(df_list, ignore_index=True)
        ## filter rows with final_score > 0.5
        combined_df = combined_df[combined_df['final_score'] > 0.9]
        combined_df = combined_df[combined_df['specificity'] < 0.01]
        print (combined_df)
        combined_df.to_csv(output_file, index=False)
        print(f"Combined host summary written to {output_file}")

if __name__ == "__main__":
    results_dir = "/home/shuaiw/borg/paper/gg_run2/"
    output_file = "./soil_all_hosts.summary.csv"
    summarize_hosts(results_dir, output_file)
