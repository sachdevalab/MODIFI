import pandas as pd

def extract_ipd_ratio(file_path):
    ## open it using pandas
    IPD_ratio_list, reverse_IPD_ratio_list = [], []
    IPD_list, reverse_IPD_list = [], []
    control_list, reverse_control_list = [], []
    df = pd.read_csv(file_path)
    # print ("csv loaded")
    ### only retain the strand == 0
    # df = df[df['strand'] == 1]
    # ### convert the modelPrediction column to list
    # # control_list = df['modelPrediction'].tolist()
    # control_list = df['tMean'].tolist()
    dict = {}
    for index, row in df.iterrows():
        if row['strand'] == 1:
            control_list.append(row['tMean'])
            if row['tpl'] not in dict:
                dict[row['tpl']] = row['tMean']
    for tpl in range(1, max(dict.keys())+1):
        if tpl not in dict:
            print (tpl)
    return control_list


file_path = "borg/split_bam_dir/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1090_C.csv"
control_list = extract_ipd_ratio(file_path)