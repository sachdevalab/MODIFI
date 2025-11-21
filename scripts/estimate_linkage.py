import numpy as np
from math import log
import pandas as pd
import os
import argparse
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import random
import math
from typing import List, Dict
# from statsmodels.stats.multitest import multipletests

from derep_motifs import MotifFilter, uniq_similar_motifs
from drep_motifs2 import motif_cluster_worker
from get_kmer_freq import kmer_freq_sim_bin_worker
# from linkage_model import compute_p_same_room
P_CUTOFF = 0.01

def specificity_weight(prevalence: float, s_max: float = 5.0) -> float:
    """
    Specificity (IDF-like): phi = clip(-ln p, 0..s_max).
    Pass prevalence already smoothed/shrunk if you have it.
    """
    # p = max(min(prevalence, 1 - 1e-9), 1e-9)
    
    p = prevalence * 0.01
    # Prevent math domain error by ensuring p > 0
    p = max(p, 1e-9)
    # print (f"prevalence: {prevalence} {p}")
    phi = -math.log(p)
    return min(max(phi, 0.0), s_max)

def linkage_score_from_counts2(motif_data, min_frac, max_sites=5000):
    """
    Adds confidence weighting based on motif site counts.
    """
    scores = []
    weights = []
    total_sites = 0

    mge_methy_sites = 0
    total_confidence = 0
    total_p_meth = 0
    probability = 1
    match_p = 1
    motif_count = 0  # Track number of motifs used in probability calculation
    miss_penalty = 0
    for m in motif_data:
        h_total = m['host_total']
        h_meth = m['host_meth']
        p_total = m['plasmid_total']
        p_meth = m['plasmid_meth']

        if h_total == 0:
            continue
        if p_total == 0:  #skip if plasmid has no motif string
            continue
        if p_total < 2: # 5
            continue

        f_host = h_meth / h_total
        f_plasmid = p_meth / p_total
        
        weight = specificity_weight(float(m["occurrence_ratio"]), s_max=6)

        if f_plasmid >= min_frac:
            motif_score = 1
            probability *= float(m["occurrence_ratio"])*0.01
            match_p *= float(m["occurrence_ratio"])*0.01
            total_sites += h_total + p_total
            motif_count += 1  # Count this motif
            mge_methy_sites += p_meth
        else:
            motif_score = 0
            probability *= (1 - float(m["occurrence_ratio"])*0.01)
            miss_penalty += 2 ** (f_host - f_plasmid) - 1
        
        scores.append(motif_score * weight)
        weights.append(weight)
        total_confidence += 10000000/m['occurrence_len'] 
        total_p_meth += p_meth


    motif_confidence = log(1+motif_count)/log(1+3)   
    motif_confidence = min(motif_confidence, 1)

    # motif_confidence = min(1, total_p_meth / 5)  

    # Confidence scaling (logarithmic)
    confidence = log(1 + total_sites) / log(1 + max_sites)
    confidence = min(confidence, 1)
    
    # final_score = linkage_score * confidence * motif_confidence
    # Normalize probability by number of motifs (geometric mean)
    if motif_count > 0:
        normalized_probability = probability ** (1.0 / motif_count)
        # normalized_probability = probability 
    else:
        normalized_probability = 1.0
    linkage_score = round(match_p, 4)
    methy_sites_confidence = min(log(1 + mge_methy_sites) / log(1 + 200), 1)
    final_score = (1 - normalized_probability) * confidence * motif_confidence 
    final_score = final_score - miss_penalty 
    final_score = max(final_score, 0)
    # final_score = -np.log(1e-16 + linkage_score)

    return {
        'specificity': round(linkage_score, 4),
        'confidence': round(confidence, 4),
        'final_score': round(final_score, 4),
        'total_sites': mge_methy_sites,
        'motif_confidence': round(motif_confidence, 4),
        'host_motif_num': len(motif_data),
    }


def linkage_score_from_counts1(motif_data, max_sites=5000):
    """
    Adds confidence weighting based on motif site counts.
    """
    scores = []
    weights = []
    total_sites = 0
    total_host_sites = 0
    restriction_signal = 1
    valid_motif_num = 0
    valid_motif_list = []
    for m in motif_data:
        h_total = m['host_total']
        h_meth = m['host_meth']
        p_total = m['plasmid_total']
        p_meth = m['plasmid_meth']

        if h_total == 0:
            continue
        if p_total == 0:  #skip if plasmid has no motif string
            continue
        if p_meth > 0:
            valid_motif_list.append(m)
            valid_motif_num += 1
        weight = h_total
        weight = min(weight, 1000)   ### set weight to 1000, once the h_total is larger than 1000
        total_sites += h_total + p_total
        total_host_sites += h_total

        f_host = h_meth / h_total

        f_plasmid = p_meth / p_total
        # motif_score = 1 - abs(f_host - f_plasmid)
        ## whether f_host and f_plasmid both are larger than 0.5
        # if f_plasmid > min_frac:
        #     motif_score = 1
        if f_plasmid > f_host:
            motif_score = 1
        else:
            motif_score = 1 - abs(f_host - f_plasmid)/f_host   ## if the f_host is only 0.5, so divided by f_host to normalize the score

        scores.append(motif_score * log(weight))
        weights.append(log(weight))
        
    if not scores:
        return {'linkage_score': 0.0, 'confidence': 0.0, 'final_score': 0.0}

    linkage_score = sum(scores) / sum(weights)

    # Confidence scaling (logarithmic)
    confidence = log(1 + total_sites) / log(1 + max_sites)
    if confidence > 1:
        confidence = 1

    # motif confidence
    # valid_motif_num = count_uniq_motif(valid_motif_list)  ## remove redundant motifs in counting
    filtered_motifs = uniq_similar_motifs(valid_motif_list)
    valid_motif_num = len(filtered_motifs)
    
    motif_confidence = log(1+valid_motif_num)/log(1+3)   
    if motif_confidence > 1:
        motif_confidence = 1
    final_score = linkage_score * confidence * motif_confidence * restriction_signal

    return {
        'linkage_score': round(linkage_score, 4),
        'confidence': round(confidence, 4),
        'final_score': round(final_score, 4),
        'total_sites': total_sites,
        'motif_confidence': round(motif_confidence, 4),
        'host_motif_num': len(motif_data)
    }

def extract_motif_data(host_df, plasmid_df, motif_cluster_dict, min_frac = 0.5, min_detect = 100):
    # Early exit if either DataFrame is empty
    if host_df.empty or plasmid_df.empty:
        return []
    
    # Use merge with inner join (default) for better performance
    motif_data = pd.merge(host_df, plasmid_df, on=['motifString', 'centerPos'], 
                         suffixes=('_host', '_plasmid'), how='inner')
    
    # Early exit if no common motifs
    if motif_data.empty:
        return []
    # Combine filtering operations into single query for better performance
    motif_data = motif_data.query(
        f'motif_modified_ratio_host > {min_frac} and motif_modified_num_host > {min_detect}'
    )
    # Early exit if no motifs pass filtering
    if motif_data.empty:
        return []
    
    # Select columns and rename in one operation
    motif_data = motif_data[['motifString', 'centerPos', 'motif_loci_num_host', 
                           'motif_modified_num_host', 'motif_loci_num_plasmid', 
                           'motif_modified_num_plasmid']]
    motif_data.columns = ['motif', 'centerPos', 'host_total', 'host_meth', 'plasmid_total', 'plasmid_meth']
    ## add a column which represents the occurrence_ratio of the motif in the host, recorded in the motif_cluster_dict
    motif_data['occurrence_ratio'] = motif_data['motif'].apply(lambda x: motif_cluster_dict[x][2])
    motif_data['cluster_id'] = motif_data['motif'].apply(lambda x: motif_cluster_dict[x][0] if x in motif_cluster_dict else None)
    motif_data['occurrence_len'] = motif_data['motif'].apply(lambda x: motif_cluster_dict[x][1] if x in motif_cluster_dict else None)

    # Convert to records - this is the most efficient way
    return motif_data.to_dict('records')

def merge_bin_profile(ctg_profile_dict, bin_name, bin_ctg_dict):
    """
    Merge the profile of all contigs in a bin into a single dataframe.
    """
    ctg_list = bin_ctg_dict[bin_name]
    ## remove the ctg not in ctg_profile_dict
    ctg_list = [ctg for ctg in ctg_list if ctg in ctg_profile_dict]
    if len(ctg_list) == 0:
        return pd.DataFrame()
    elif len(ctg_list) == 1:
        return ctg_profile_dict[ctg_list[0]]
    else:
        bin_df = ctg_profile_dict[ctg_list[0]]
        ## only retain motifString,centerPos,motif_loci_num, motif_modified_num, and motif_modified_ratio
        bin_df = bin_df[['motifString', 'centerPos', 'motif_loci_num', 'motif_modified_num', 'motif_modified_ratio']]
        for ctg_name in ctg_list[1:]:
            if ctg_name not in ctg_profile_dict:
                print (f"{ctg_name} has no methylation profile.")
                continue
            df = ctg_profile_dict[ctg_name]
            merge_df = pd.merge(bin_df, df, on = ['motifString', 'centerPos'], suffixes=('_x', '_y'))
            merge_df['motif_loci_num'] = merge_df['motif_loci_num_x'] + merge_df['motif_loci_num_y']
            merge_df['motif_modified_num'] = merge_df['motif_modified_num_x'] + merge_df['motif_modified_num_y']
            ## avoid divided by 0, keep original value
            merge_df['motif_modified_ratio'] = merge_df['motif_modified_num'] / merge_df['motif_loci_num'].replace(0, 1)
            bin_df = merge_df[['motifString', 'centerPos', 'motif_loci_num', 'motif_modified_num', 'motif_modified_ratio']]
        return bin_df

def merge_bin_motif_file(ctg_motif_dict, bin_name, bin_ctg_dict, min_frac, min_detect):
    ctg_list = bin_ctg_dict[bin_name]
    if len(ctg_list) == 1:
        if ctg_list[0] not in ctg_motif_dict:
            return pd.DataFrame()
        return ctg_motif_dict[ctg_list[0]]
    else:
        motif_df = pd.DataFrame()
        for ctg_name in ctg_list:
            if ctg_name not in ctg_motif_dict:
                # print (f"{ctg_name} has no motif file.")
                continue
            df = ctg_motif_dict[ctg_name]
            if not df.empty:
                motif_df = pd.concat([motif_df, df], axis=0)
        
        # If no motif data found, return empty DataFrame
        if motif_df.empty:
            return motif_df
            
        aggregated_motif = motif_df.groupby(['motifString', 'centerPos']).agg({
            'nGenome': 'sum',
            'nDetected': 'sum'
        }).reset_index()
        
        # Recalculate the modified ratio
        aggregated_motif['fraction'] = (
            aggregated_motif['nDetected'] / 
            aggregated_motif['nGenome'].replace(0, 1)
        )
        ## only retain the motifs with fraction >= min_frac and nDetected >= MIN_DETECT
        aggregated_motif = aggregated_motif[
            (aggregated_motif['fraction'] > min_frac) & 
            (aggregated_motif['nDetected'] > min_detect)
        ]
        
        return aggregated_motif

def drep_bin_motif(bin_motif, motif_cluster_dict):
    ### if two motif has the same cluster ID, then only keep one with highest nDetected
    bin_motif['cluster_id'] = bin_motif.apply(
        lambda row: motif_cluster_dict.get(row['motifString'])[0] if row['motifString'] in motif_cluster_dict else None, axis=1
    )
    filtered_motifs = []
    seen_clusters = set()
    for _, row in bin_motif.sort_values(by='nDetected', ascending=False).iterrows():
        cluster_id = row['cluster_id']
        if cluster_id and cluster_id not in seen_clusters:
            filtered_motifs.append(row)
            seen_clusters.add(cluster_id)
        elif not cluster_id:  # Keep motifs without a cluster ID
            filtered_motifs.append(row)
    return pd.DataFrame(filtered_motifs)


def merge_bin_motif(bin_cov_dict, bin_ctg_dict, ctg_motif_dict, ctg_profile_dict, motif_cluster_dict, min_frac, min_detect):
    bin_df_dict = {}
    bin_motif_dict = {}
    for bin_name in bin_cov_dict:
        bin_df = merge_bin_profile(ctg_profile_dict, bin_name, bin_ctg_dict)
        bin_motif = merge_bin_motif_file(ctg_motif_dict, bin_name, bin_ctg_dict, min_frac, min_detect)
        bin_motif = drep_bin_motif(bin_motif, motif_cluster_dict)
        if len(bin_df) == 0 and len(bin_motif) == 0:
            print(f"Bin {bin_name} has no methylation profile or motif data.")
            continue
        bin_df_dict[bin_name] = bin_df
        bin_motif_dict[bin_name] = bin_motif
    return bin_df_dict, bin_motif_dict
        
def estimate_cov(cov_dict, bin_name, bin_ctg_dict):
    ctg_cov_list = []
    for ctg in bin_ctg_dict[bin_name]:
        if ctg in cov_dict:
            ctg_cov_list.append(cov_dict[ctg])
    if len(ctg_cov_list) >0:
        bin_cov = round(np.mean(ctg_cov_list),2)
    else:
        bin_cov = 'NA'
    return bin_cov

def summary_motif_info(motif_data):
    ## summary it into a string, use ; separate motifs, use : to separate the motif info
    motif_info = []
    for m in motif_data:
        motif_info.append(m['motif'] + ":" + str(m['centerPos']) + ":" + str(m['host_total']) + \
                           ":" + str(m['host_meth']) + ":" + str(m['plasmid_total']) + ":" + \
                            str(m['plasmid_meth']) + ":" + str(m['occurrence_ratio']) + ":" + str(m['occurrence_len']))
    motif_info = ";".join(motif_info)
    return motif_info

def bin_worker(bin_df, plasmid_df, bin_motif, min_frac, bin_name, min_detect, motif_cluster_dict):
    motif_data = extract_motif_data(bin_df, plasmid_df, motif_cluster_dict, min_frac, min_detect)
    motif_data = filter_motifs(bin_motif, motif_data)
    # motif_filter = MotifFilter(motif_data)
    # motif_data = motif_filter.filter()

    result = linkage_score_from_counts2(motif_data, min_frac)
    return result, motif_data, bin_name

def for_each_plasmid(bin_df_dict, bin_motif_dict, plasmid_name, profile_dir, host_dir,
                      cov_dict, bin_cov_dict, motif_cluster_dict, bin_ctg_dict, threads, min_frac = 0.5, min_detect = 100):
    import gc
    plasmid_profile = f"{profile_dir}/{plasmid_name}.motifs.profile.csv"
    # cov_dict = load_coverage(host_dir)
    if plasmid_name not in cov_dict:
        MGE_cov = 'NA'
    else:
        MGE_cov = cov_dict[plasmid_name]
    ## check if the plasmid profile exists
    if not os.path.exists(plasmid_profile):
        print (f"plasmid_profile {plasmid_profile} does not exist.")
        return []
    plasmid_df = pd.read_csv(plasmid_profile)

    # # Get all valid bins to process
    # valid_bins = []
    # for bin_name in bin_cov_dict:
    #     bin_df = bin_df_dict[bin_name]
    #     if len(bin_df) > 0:
    #         valid_bins.append(bin_name)
    valid_bins = list(bin_df_dict.keys())
    total_bins = len(valid_bins)
    if total_bins == 0:
        print(f"No valid bins found for plasmid {plasmid_name}")
        return []
    
    print(f"Processing {total_bins} bins for plasmid {plasmid_name}")
    
    all_result = []
    batch_size = 1000  # Process bins in batches
    effective_threads = min(threads, 64)  # Limit threads for bin processing
    
    for i in range(0, total_bins, batch_size):
        batch_bins = valid_bins[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (total_bins + batch_size - 1) // batch_size
        
        print(f"Processing bin batch {batch_num}/{total_batches} ({len(batch_bins)} bins) for plasmid {plasmid_name}")

        with ProcessPoolExecutor(max_workers=effective_threads) as executor:
            batch_futures = []
            for bin_name in batch_bins:
                bin_df = bin_df_dict[bin_name]
                bin_motif = bin_motif_dict[bin_name]

                future = executor.submit(
                    bin_worker,
                    bin_df = bin_df,
                    plasmid_df = plasmid_df,
                    bin_motif = bin_motif,
                    min_frac = min_frac,
                    bin_name = bin_name,
                    min_detect = min_detect,
                    motif_cluster_dict = motif_cluster_dict,
                )
                batch_futures.append(future)

            # Process batch results
            batch_results = []
            for future in tqdm(as_completed(batch_futures), 
                             total=len(batch_futures), 
                             desc=f"Plasmid {plasmid_name} Batch {batch_num}"):
                try:
                    result = future.result(timeout=60)  # 60-second timeout per bin
                    batch_results.append(result)
                except Exception as e:
                    print(f"Error processing bin in batch {batch_num}: {e}")
        
        # Add batch results to main results
        all_result.extend(batch_results)
        
        # Force garbage collection between batches
        gc.collect()
        
        print(f"Bin batch {batch_num} completed for plasmid {plasmid_name}. Total results so far: {len(all_result)}")

    data = []
    motif_data_dict = {}
    i = 0
    for each_result in all_result:
        result, motif_data, bin_name = each_result
        if not motif_data:
            continue
        result['host'] = bin_name
        result['host_cov'] = bin_cov_dict[bin_name]
        result['MGE_cov'] = MGE_cov
        result['MGE'] = plasmid_name
        result['motif_info'] = summary_motif_info(motif_data)
        motif_data_dict[result['host']] = motif_data
        data.append(result)
        i += 1

    print (f"### Processed bins for {plasmid_name}, total {len(data)} results.")
    ## convert the data to a df, and sort by final_score
    data = pd.DataFrame(data)
    final_score_list = []
    if len(data) > 0:
        ## calculate a pvalue to indicate proportion of the final_score that are larger than the value for this plasmid only
        data['self_pvalue']  = data['final_score'].apply(lambda x: sum(data['final_score'] >= x) / len(data) )
        # Separate rows with specificity < P_CUTOFF and others
        low_specificity = data[data['specificity'] < P_CUTOFF].sort_values(by='final_score', ascending=False, ignore_index=True)
        high_specificity = data[data['specificity'] >= P_CUTOFF].sort_values(by='final_score', ascending=False, ignore_index=True)

        # Combine them: low specificity first, then high specificity
        data = pd.concat([low_specificity, high_specificity], ignore_index=True)
        
        ## if more than one hosts have ths highest final_score, then sort them by cos_sim, and keep them all
        data = sort_top_by_cos_sim(data, bin_ctg_dict, os.path.join(host_dir, '../'))
        ## host column the first, final_score the second, specificity the third, confidence the forth, total_sites the fifth
        data = data[['MGE', 'host', 'final_score', 'specificity', 'self_pvalue', 'confidence', 'motif_confidence', 'total_sites', 'host_motif_num', 'MGE_cov', 'host_cov','motif_info']]
        data['self_pvalue'] = data['self_pvalue'].round(4)
        
        final_score_list = data['final_score'].tolist()
        print (data.head(5))
        ## output the data to a csv file
        host_prediction = os.path.join(host_dir, f"{plasmid_name}.host_prediction.csv")
        data.to_csv(host_prediction, index = False)
    ## output the motif data to a csv file
    motif_data_file = os.path.join(host_dir, f"{plasmid_name}.motif_data.csv")
    with open(motif_data_file, 'w') as f:
        for index, rows in data.iterrows():
            print (rows.to_dict(), motif_data_dict[rows['host']], file = f)
    return final_score_list

def sort_top_by_cos_sim(data, bin_ctg_dict, work_dir):
    ## check if more than one hosts have ths highest final_score, then sort them by cos_sim, and keep them all
    top_final_score = data.iloc[0]['final_score']
    top_data = data[data['final_score'] == top_final_score].copy()
    other_data = data[data['final_score'] < top_final_score].copy()

    top_data['cos_sim'] = 0.0
    if len(top_data) > 1 and len(top_data) < 4 and top_final_score > 0:
        ## calculate cos_sim for these hosts
        for i, row in top_data.iterrows():
            if row['MGE'] not in bin_ctg_dict:
                bin_ctg_dict[row['MGE']] = [row['MGE']]
            bin_name1, bin_1_gc, bin_2_gc, cos_sim = kmer_freq_sim_bin_worker(
                row['MGE'], row['host'], bin_ctg_dict, work_dir
            )
            top_data.loc[i, 'cos_sim'] = cos_sim
        ## sort the top_data by cos_sim
        top_data = top_data.sort_values(by = 'cos_sim', ascending = False, ignore_index = True)
        return pd.concat([top_data, other_data], ignore_index=True)
    else:
        return data

def report_gc(data, host_dir, bin_ctg_dict, threads):
    """
    report GC for MGE and host
    also report the cosine similarity for MGE and host tetranuceleotide frequency
    """
    print (f"report_gc for {len(data)} MGE-host pairs.")
    work_dir = os.path.join(host_dir, "../")
    ## add new columns for MGE_gc, host_gc, cos_sim
    data['MGE_gc'] = 0.0
    data['host_gc'] = 0.0
    data['cos_sim'] = 0.0

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = []
        for index, row in data.iterrows():
            if row['MGE'] not in bin_ctg_dict:
                bin_ctg_dict[row['MGE']] = [row['MGE']]
            future = executor.submit(
                kmer_freq_sim_bin_worker,
                row['MGE'], row['host'], bin_ctg_dict, work_dir
            )
            futures.append(future)

    results = []
    for future in tqdm(as_completed(futures), total=len(futures)):
        result = future.result()
        if result is not None:
            results.append(result)
    for idx, result in enumerate(results):
        MGE, MGE_gc, host_gc, tetra_sim = result
        # print(MGE_gc, host_gc, tetra_sim)
        # Use MGE as index to update the DataFrame
        if MGE in data['MGE'].values:
            data.loc[data['MGE'] == MGE, 'MGE_gc'] = MGE_gc
            data.loc[data['MGE'] == MGE, 'host_gc'] = host_gc
            data.loc[data['MGE'] == MGE, 'cos_sim'] = tetra_sim
    # data = data[['MGE', 'host', 'final_score', 'linkage_score', 'confidence', 'motif_confidence', 'total_sites',\
    #               'host_motif_num', 'MGE_gc', 'host_gc', 'cos_sim', 'MGE_cov', 'host_cov','motif_info']]
    return data

def read_genomad(genomad_file):
    MGE_dict = {}
    print (f"Reading {genomad_file}...")
    genomad = pd.read_csv(genomad_file, sep = "\t")
    for i, row in genomad.iterrows():
        
        if re.search('\|provirus', row['seq_name']):
            continue
        if row['seq_name'] == 'seq_name':
            continue
        # print (row['seq_name'])
        MGE_dict[row['seq_name']] = row
    return MGE_dict

def filter_motifs(host_motif, motif_data):
    # Early exit if no data
    if host_motif.empty or not motif_data:
        return []
    
    # Use vectorized operations instead of iterrows() - much faster
    # Create tags using vectorized string operations
    host_motif_tags = (host_motif['motifString'] + "_tga_" + host_motif['centerPos'].astype(str)).tolist()
    motif_tag_set = set(host_motif_tags)  # Use set for O(1) lookup instead of dict
    
    # Filter motifs using list comprehension - faster than loop + append
    new_motif_data = [
        m for m in motif_data 
        if (m['motif'] + "_tga_" + str(m['centerPos'])) in motif_tag_set
    ]
    
    return new_motif_data

def count_MGE_with_motif(plasmid_name, profile_dir):
    plasmid_motif_file = os.path.join(profile_dir,"../motifs", f"{plasmid_name}.motifs.csv")
    if not os.path.exists(plasmid_motif_file):
        print (f"plasmid_motif_file {plasmid_motif_file} does not exist.")
        return 0
    plasmid_motif = pd.read_csv(plasmid_motif_file)
    return len(plasmid_motif)

def load_coverage(host_dir):
    cov_file = os.path.join(host_dir, "../mean_depth.csv")
    if not os.path.exists(cov_file):
        print (f"{cov_file} does not exist.")
        return {}
    cov_df = pd.read_csv(cov_file)
    cov_dict = {}
    for i, row in cov_df.iterrows():
        contig = row['contig']
        mean_depth = row['depth']
        cov_dict[contig] = mean_depth
    return cov_dict

def summary_host(host_dir, bin_ctg_dict, threads, all_final_score_list, MGE_dict, n_iter = 10000):
    data = []
    for file in os.listdir(host_dir):
        if file.endswith(".host_prediction.csv"):
            # plasmid_name = file.split(".")[0]
            host_prediction = os.path.join(host_dir, file)
            ## check if host_prediction is empty
            if os.path.getsize(host_prediction) == 0:
                continue
            df = pd.read_csv(host_prediction)
            ## extract the first row
            ## add new column for plasmid_name at the start
            if len(df) > 0:
                # df.insert(0, 'MGE', plasmid_name)
                best_host = df.iloc[0].copy()  # Make a copy to avoid SettingWithCopyWarning
                if best_host['MGE'] not in MGE_dict:
                    # print (f"{best_host['MGE']} not in MGE_dict, skip.")
                    continue

                if best_host['MGE'] in MGE_dict and 'length' in MGE_dict[best_host['MGE']]:
                    best_host['MGE_len'] = MGE_dict[best_host['MGE']]['length']
                else:
                    best_host['MGE_len'] = 0
                data.append(best_host)  # Append the modified copy

    ## randomly select n_iter values from all_final_score_list
    print (f"{len(all_final_score_list)} final scores in all_final_score_list used as background.")
    num = int(min(n_iter, 1 *  len(all_final_score_list)))
    selected_scores = random.sample(all_final_score_list, num)
    # print (sorted(selected_scores, reverse=True))

    data = pd.DataFrame(data)

    ## sort by final_score
    if len(data) > 0:
        # data = data.sort_values(by = 'final_score', ascending = False, ignore_index = True)
        low_specificity = data[data['specificity'] < P_CUTOFF].sort_values(by='final_score', ascending=False, ignore_index=True)
        high_specificity = data[data['specificity'] >= P_CUTOFF].sort_values(by='final_score', ascending=False, ignore_index=True)
        
        # Combine them: low specificity first, then high specificity
        data = pd.concat([low_specificity, high_specificity], ignore_index=True)

        data = report_gc(data, host_dir, bin_ctg_dict, threads)
        ## add pvalue for final_score
        data['pvalue'] = data['final_score'].apply(lambda x: sum(1 for score in selected_scores if score >= x) / len(selected_scores))
        data['pvalue'] = data['pvalue'].round(4)
        # print (data["MGE_gc"])
        ## resort the columns
        data = data[['MGE', 'MGE_len', 'host', 'final_score', 'specificity', 'pvalue', 'self_pvalue', 'MGE_gc', 'host_gc', 'cos_sim', 
                    'MGE_cov', 'host_cov', 'host_motif_num', 'confidence', 
                    'motif_confidence',  'total_sites', 'motif_info']]
    ## output the data to a csv file
    host_summary = os.path.join(host_dir, "../", "host_summary.csv")
    data.to_csv(host_summary, index = False)
    filter_linkage(data, host_dir)

def filter_linkage(data, host_dir, small_mge_dp = 10):
    filter_host_summary = os.path.join(host_dir, "../", "host_summary.filter.csv")
    if len(data) > 0:
        data = data[data['pvalue'] < 0.05]
        data = data[data['final_score'] > 0.5]
    ## if MGE_len < 5000, the min depth should be 10 for both host and MGE
    filtered_data = []
    for index, row in data.iterrows():
        if row['MGE_len'] < 5000:
            if row['MGE_cov'] < small_mge_dp or row['host_cov'] < small_mge_dp:
                continue
        filtered_data.append(row)
    filtered_data = pd.DataFrame(filtered_data)
    filtered_data.to_csv(filter_host_summary, index = False)

def get_bin_cov(cov_dict, bin_ctg_dict, min_ctg_cov, whole_ref, MGE_dict):
    whole_ref_fai = whole_ref + ".fai"
    if not os.path.exists(whole_ref_fai):
        raise FileNotFoundError(f"{whole_ref_fai} does not exist. Please provide a valid whole reference fasta file.")
    ctg_len_dict = {}
    with open(whole_ref_fai, 'r') as fai_file:
        for line in fai_file:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                ctg_name = parts[0]
                ctg_len = int(parts[1])
                ctg_len_dict[ctg_name] = ctg_len

    bin_cov_dict = {}
    for bin_name in bin_ctg_dict:
        if bin_name in MGE_dict:
            continue
        ctg_cov_list = []
        bin_length = 0
        for ctg in bin_ctg_dict[bin_name]:
            bin_length += ctg_len_dict[ctg]
            if ctg in cov_dict:
                ctg_cov_list.append(cov_dict[ctg] * ctg_len_dict[ctg])

        if len(ctg_cov_list) >0:
            # bin_cov = round(np.mean(ctg_cov_list),2)
            bin_cov = round(np.sum(ctg_cov_list)/bin_length,2)
            if bin_cov >= min_ctg_cov:
                bin_cov_dict[bin_name] = bin_cov
        # else:
        #     bin_cov = 'NA'
    return bin_cov_dict

def batch_MGE_invade(plasmid_file, profile_dir, host_dir, whole_ref, bin_file=None, min_frac = 0.5, threads = 1, min_ctg_cov = 5, min_detect = 100):
    ### cluster motifs
    fai = whole_ref + ".fai"
    if not os.path.exists(fai):
        raise FileNotFoundError(f"{fai} does not exist. Please provide a valid whole reference fasta file.")
    output_dir = os.path.join(host_dir, "../")
    motif_file = os.path.join(output_dir, "motif_profile.csv")
    length_df = motif_cluster_worker(motif_file, fai, output_dir, min_frac)
    ## get a dict to record the motif cluster of each motif
    motif_cluster_dict = {}
    for i, row in length_df.iterrows():
        motif = row['motifString']
        cluster = row['cluster_id']
        occurrence_length = row['total_length']
        occurrence_ratio = row['percentage_of_profile']
        motif_cluster_dict[motif] = [cluster, occurrence_length, occurrence_ratio]

    MGE_dict = read_genomad(plasmid_file)
    cov_dict = load_coverage(host_dir)
    ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict = load_ctg_motifs_parallele(profile_dir, threads)
    # ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict = load_ctg_motifs_single(profile_dir, MGE_dict)
    if bin_file is not None:
        if not os.path.exists(bin_file):
            print (f"{bin_file} does not exist.")
        ctg_bin_dict, bin_ctg_dict = load_bin(bin_file)
    else:
        ctg_bin_dict, bin_ctg_dict = ctg_single_dict, single_ctg_dict
    print (f"Loading is done, {len(bin_ctg_dict)} bins.")
    # """
    bin_cov_dict = get_bin_cov(cov_dict, bin_ctg_dict, min_ctg_cov, whole_ref, MGE_dict)  ## estimate the coverage for each bin, and remove bins with cov < min_ctg_cov
    print (f"Estimated {len(bin_cov_dict)} bins with coverage >= {min_ctg_cov}.")
    
    bin_df_dict, bin_motif_dict = merge_bin_motif(bin_cov_dict, bin_ctg_dict, ctg_motif_dict, ctg_profile_dict, motif_cluster_dict, min_frac, min_detect)
    print (f"Loaded {len(bin_df_dict)} bins with methylation profile and motifs.")
    
    print ("threads number for linkage prediction: ", threads)
    ## remove the MGE with cov < min_ctg_cov
    MGE_dict = {k: v for k, v in MGE_dict.items() if k in cov_dict and cov_dict[k] >= min_ctg_cov}
    print ("MGE number: ", len(MGE_dict))

    i = 0
    all_final_score_list = []
    for plasmid_name in MGE_dict:
        # MGE_motif_num = count_MGE_with_motif(plasmid_name, profile_dir)
        # if MGE_motif_num == 0:  # to speed up
        #     print (f"Skip {plasmid_name} with {MGE_motif_num} motifs.")
        #     continue
        # print (f"Processing {i}-th/{len(MGE_dict)} {plasmid_name} with {MGE_motif_num} original motifs.")
        # if plasmid_name != "ocean_1_1355_L":
        #     continue
        print (f"Processing {i}-th/{len(MGE_dict)} {plasmid_name}.")

        final_score_list = for_each_plasmid(
            bin_df_dict = bin_df_dict,
            bin_motif_dict = bin_motif_dict,
            plasmid_name = plasmid_name,
            profile_dir = profile_dir,
            host_dir = host_dir,
            cov_dict = cov_dict,
            bin_cov_dict = bin_cov_dict,
            threads = threads,
            min_frac = min_frac,
            min_detect = min_detect,
            motif_cluster_dict = motif_cluster_dict,
            bin_ctg_dict = bin_ctg_dict,
        )
        i += 1
        if final_score_list is not None:
            all_final_score_list += final_score_list


    summary_host(host_dir, bin_ctg_dict, threads, all_final_score_list,MGE_dict)

def load_bin(bin_file):
    bin_df = pd.read_csv(bin_file, sep= "\t", header = 0)
    ## add header for the bin_df
    bin_df.columns = ['contig', 'bin_id']
    # print (bin_df)
    ctg_bin_dict = {}
    bin_ctg_dict = defaultdict(list)

    for i, row in bin_df.iterrows():
        contig = row['contig']
        bin_id = row['bin_id']
        ctg_bin_dict[contig] = bin_id
        bin_ctg_dict[bin_id].append(contig)
    return ctg_bin_dict, bin_ctg_dict

def process_file_with_args(args):
    file, profile_dir = args
    try:
        if file.endswith(".motifs.profile.csv"):
            ctg_name = file.split(".motifs.profile.csv")[0]
            host_motif = os.path.join(profile_dir, "../motifs", f"{ctg_name}.motifs.csv")
            if not os.path.exists(host_motif):
                print(f"process_file_with_args host_motif {host_motif} does not exist.")
                return None
            
            host_profile = os.path.join(profile_dir, file)
            if not os.path.exists(host_profile):
                print(f"host_profile {host_profile} does not exist.")
                return None
            
            # Use faster CSV reading with engine='c'
            motif_df = pd.read_csv(host_motif, engine='c')
            host_df = pd.read_csv(host_profile, engine='c')
            
            return {
                "ctg_name": ctg_name,
                "motif_df": motif_df,
                "host_df": host_df
            }
    except Exception as e:
        print(f"Error processing file {file}: {e}")
        return None
    return None

def load_ctg_motifs_parallele(profile_dir, threads=4):
    import gc
    print(f"load_ctg_motifs from {profile_dir}...")
    
    # Get all files to process
    all_files = [file for file in os.listdir(profile_dir) if file.endswith(".motifs.profile.csv")]
    # all_files = ["cow_1_3388_L.motifs.profile.csv"]
    total_files = len(all_files)
    
    if total_files == 0:
        print("No motif profile files found.")
        return {}, defaultdict(list), {}, {}
    
    print(f"Found {total_files} files to process")
    
    # Initialize result containers
    results = []
    processed_count = 0
    
    # Process files in batches to avoid memory issues
    batch_size = 100
    
    # Use fewer workers for better memory management
    effective_threads = min(threads, 64)
    
    for i in range(0, total_files, batch_size):
        batch_files = all_files[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (total_files + batch_size - 1) // batch_size
        
        print(f"Processing batch {batch_num}/{total_batches} ({len(batch_files)} files)")
        
        # Process current batch
        with ProcessPoolExecutor(max_workers=effective_threads) as executor:
            # Create futures for current batch only
            batch_futures = []
            for file in batch_files:
                future = executor.submit(
                    process_file_with_args,
                    (file, profile_dir)
                )
                batch_futures.append(future)
            
            # Process results as they complete
            batch_results = []
            for future in tqdm(as_completed(batch_futures), 
                             total=len(batch_futures), 
                             desc=f"Batch {batch_num}"):
                try:
                    result = future.result(timeout=30)  # 30-second timeout per file
                    if result is not None:
                        if result["motif_df"].empty or result["host_df"].empty:
                            # print(f"Skipping empty result for {result['ctg_name']}")
                            continue
                        batch_results.append(result)
                except Exception as e:
                    print(f"Error processing file in batch {batch_num}: {e}")
                
                processed_count += 1
                
                # Log progress every 50 files
                if processed_count % 50 == 0:
                    print(f"Processed {processed_count}/{total_files} files, "
                          f"found {len(results) + len(batch_results)} valid results")
        
        # Add batch results to main results
        results.extend(batch_results)
        
        # Force garbage collection between batches to free memory
        gc.collect()
        
        print(f"Batch {batch_num} completed. Total valid results so far: {len(results)}")

    print(f"All batches completed. Processing {len(results)} valid results into dictionaries...")

    # Combine results into dictionaries
    ctg_single_dict = {}
    single_ctg_dict = defaultdict(list)
    ctg_profile_dict = {}
    ctg_motif_dict = {}

    for result in results:
        ctg_name = result["ctg_name"]
        ctg_single_dict[ctg_name] = ctg_name
        single_ctg_dict[ctg_name].append(ctg_name)
        ctg_profile_dict[ctg_name] = result["host_df"]
        ctg_motif_dict[ctg_name] = result["motif_df"]
    
    print(f"load_ctg_motifs done. Processed {len(results)} contigs successfully.")
    return ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get invasion score of MGE", add_help=False, \
    usage="%(prog)s -h", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("--work_dir", type=str, help="<str> work_dir", metavar="\b")
    # required.add_argument("--all_profiles", type=str, help="<str> separate by space", nargs="+", metavar="\b")
    optional.add_argument("--plasmid_file", type=str, help="<str> *_plasmid_summary.tsv by genomad.", metavar="\b")
    optional.add_argument("--plasmid", type=str, help="<str> plasmid contig name.", metavar="\b")
    optional.add_argument("--bin_file", type=str, help="<str> binning file: contig bin_id", metavar="\b")
    optional.add_argument("--min_frac", type=float, help="<str> output motif summary.", default=0.5, metavar="\b")
    ## add threads
    optional.add_argument("--threads", type=int, help="<int> number of threads.", default=1, metavar="\b")

    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())



    min_frac = args["min_frac"]
    # plasmid_contig = "E_coli_H10407_4"
    # work_dir = "/home/shuaiw/borg/bench/zymo6_NM3/"
    # genomad_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_genomad/merged2_summary/merged2_plasmid_summary.tsv"
    work_dir = args["work_dir"]
    profile_dir = os.path.join(work_dir, "profiles")
    host_dir = os.path.join(work_dir, "hosts")
    bin_file = args["bin_file"]
    os.makedirs(host_dir, exist_ok = True)


    if args["plasmid_file"]:
        batch_MGE_invade(args["plasmid_file"], profile_dir, host_dir, bin_file, min_frac, args["threads"])
    elif args["plasmid"]:
        plasmid_name = args["plasmid"]
        MGE_motif_num = count_MGE_with_motif(plasmid_name, profile_dir)
        print (f"Processing {plasmid_name} with {MGE_motif_num} original motifs.")
        ctg_single_dict, single_ctg_dict, ctg_profile_dict, ctg_motif_dict = load_ctg_motifs_parallele(profile_dir, {}, args["threads"])
        print (f"Loading is done, {len(single_ctg_dict)} contigs.", len(ctg_motif_dict))
        if bin_file is not None:
            if not os.path.exists(bin_file):
                print (f"bin_file {bin_file} does not exist.")
            ctg_bin_dict, bin_ctg_dict = load_bin(bin_file)
        else:
            ctg_bin_dict, bin_ctg_dict = ctg_single_dict, single_ctg_dict
        bin_df_dict, bin_motif_dict = merge_bin_motif(bin_ctg_dict, ctg_motif_dict, ctg_profile_dict)
        for_each_plasmid(bin_df_dict, bin_motif_dict, bin_ctg_dict, ctg_profile_dict, ctg_motif_dict, args["plasmid"], profile_dir, host_dir, min_frac, {})
    else:
        print ("Please provide either --plasmid_file or --plasmid.")

    # extract the info saved in the host_dir
    summary_host(host_dir)



    """

    min_frac = 0.5
    # host_profile = "/home/shuaiw/borg/bench/zymo6_NM3/profiles/E_coli_H10407_1.motifs.profile.csv"
    host_profile = "/home/shuaiw/borg/all_test_ccs3/profiles/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_4094_L.motifs.profile.csv"
    host_motif = "/home/shuaiw/borg/all_test_ccs3/motifs/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_4094_L.motifs.csv"
    plasmid_profile = "/home/shuaiw/borg/all_test_ccs3/profiles/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META_1354_L.motifs.profile.csv"

    motif_data = extract_motif_data(host_profile, plasmid_profile, min_frac)
    print (motif_data)
    motif_data = filter_motifs(host_motif, motif_data)
    print (motif_data)

    result = linkage_score_from_counts(motif_data)
    print(result)"
    """

# python cal_linkage_score.py --work_dir /home/shuaiw/borg/bench/zymo_new_ref --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2_genomad/merged2_summary/merged2_plasmid_summary.tsv
# python cal_linkage_score.py  --work_dir /home/shuaiw/borg/all_test_ccs3 --plasmid_file /home/shuaiw/borg/contigs/borg_pure.txt
# python cal_linkage_score.py  --work_dir /home/shuaiw/borg/bench/all_ccs_1k --plasmid_file /home/shuaiw/borg/contigs/borg_pure.txt
# python cal_linkage_score.py  --work_dir /home/shuaiw/borg/bench/all_ccs_1k --plasmid_file /home/shuaiw/borg/contigs/genomad/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_summary/SR-VP_9_9_2021_81_5A_0_75m_PACBIO-HIFI_HIFIASM-META.contigs_plasmid_summary.tsv
# python cal_linkage_score.py --work_dir /home/shuaiw/borg/bench/zymo_new_ref_p0.05_cov1_s30 --plasmid_file /home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list

