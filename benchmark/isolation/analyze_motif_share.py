import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform, cosine
from sklearn.metrics.pairwise import cosine_similarity
import re
import sys
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'motif_change'))
from sample_object import My_sample, Isolation_sample, get_unique_motifs,My_contig
from check_motif_change import read_drep_cluster, clean_profile

relation_order = ['same_isolate', 'same_strain','same_species', 'same_genus', 'same_family', 'same_order', 
                    'same_class', 'same_phylum', 'same_domain']

def read_motif_freq(motif_freq_file, prefix, all_motif_freq):
    motif_freq_df = pd.read_csv(motif_freq_file)
    ## make the column percentage_of_profile * 0.01
    motif_freq_df["percentage_of_profile"] = motif_freq_df["percentage_of_profile"] * 0.01
    ## only keep the rows with percentage_of_profile > 0.01
    # motif_freq_df = motif_freq_df[motif_freq_df["percentage_of_profile"] > 0.04]
    ## remove the rows with motifString's and start end letter is not in GATC
    motif_freq_df = motif_freq_df[motif_freq_df["motifString"].str.endswith(('G', 'A', 'T', 'C'))]
    motif_freq_df = motif_freq_df[motif_freq_df["motifString"].str.startswith(('G', 'A', 'T', 'C'))]
    ## percentage_of_profile > 0.1, set it as 0.1
    # motif_freq_df.loc[motif_freq_df["percentage_of_profile"] > 0.1, "percentage_of_profile"] = 0.1
    motif_freq_df["prefix"] = prefix
    for i, row in motif_freq_df.iterrows():
        all_motif_freq.append([row["percentage_of_profile"], prefix, row["motifString"]])
    return motif_freq_df, all_motif_freq

def read_motif_freq_ctg(profile, all_motif_freq, mge_dict):
    df = pd.read_csv(profile)
    ## get the column names
    cols = df.columns.tolist()
    for i, row in df.iterrows():
        # print (row)
        for j in range(1, len(cols)):
            if cols[j] in mge_dict:
                ctg_type = "MGE"
            else:
                ctg_type = "Host"
            all_motif_freq.append([row[cols[j]], cols[j], row["motifString"], ctg_type])
    return all_motif_freq

def quantify_sharing(sample_obj, mge_dict, depth_dict, length_dict, unique_motifs, represent_contig_list):

    df = pd.read_csv(sample_obj.profile)
    contig_cols = df.columns.tolist()[1:]  # Skip the motifString column
    
    contig_annotations = []
    for contig in set(contig_cols):
        # Determine contig type
        if contig in mge_dict:
            ctg_type = "MGE"
            mge_type = mge_dict[contig]
        else:
            if contig not in represent_contig_list:
                ## only keep the host contigs in represent_contig_list
                continue
            ctg_type = "Host"
            mge_type = "N/A"
        
        # Get depth and length
        depth = depth_dict.get(contig, 0.0)
        length = length_dict.get(contig, 0)
        
        contig_annotations.append({
            'contig': contig,
            'contig_type': ctg_type,
            'mge_type': mge_type,
            'depth': depth,
            'length': length,
            'prefix': sample_obj.prefix,
            'phylum': sample_obj.phylum,
            'lineage': sample_obj.lineage
        })
    
    # Convert to DataFrame
    annotations_df = pd.DataFrame(contig_annotations)
    
    # Transform the motif profile from wide to long format
    df_melted = df.melt(id_vars=['motifString'], 
                        var_name='contig', 
                        value_name='motif_frequency')
    ## add motif string number onto df_melted
    for i, row in df_melted.iterrows():
        motif_str = row['motifString']
        contig = row['contig']
        ctg_obj = My_contig(sample_obj.prefix, sample_obj.all_dir , contig, data_type="isolation")
        motif_loci_num_dict = ctg_obj.read_ctg_profile()
        # print ("####", contig, motif_str, motif_loci_num_dict)
        if motif_str in motif_loci_num_dict:
            loci_num = motif_loci_num_dict[motif_str]
            df_melted.at[i, 'motif_loci_num'] = loci_num
    # Merge with annotations
    df_enhanced = df_melted.merge(annotations_df, on='contig', how='left')
    
    df_enhanced = df_enhanced[df_enhanced['motifString'].isin(unique_motifs)]
    # filter annotations_df by sample_obj.depth_cutoff and sample_obj.length_cutoff
    annotations_df = annotations_df[
        (annotations_df['depth'] >= sample_obj.depth_cutoff) &
        (annotations_df['length'] >= sample_obj.length_cutoff)
    ]
    jaccard_scores , present_motifs, filtered_mge_num = cal_jaccard(df_enhanced, sample_obj.depth_cutoff, sample_obj.length_cutoff)
    return jaccard_scores, present_motifs, annotations_df, filtered_mge_num
    # return df_enhanced

def cal_jaccard(df_enhanced, depth_cutoff, length_cutoff, bin_freq = 0.3):
    ## filter the dataframe based on depth and length cutoff
    filtered_df = df_enhanced[
        (df_enhanced['depth'] >= depth_cutoff) &
        (df_enhanced['length'] >= length_cutoff)
    ].copy()

    ## binary the motif_frequency based on bin_freq
    filtered_df['motif_frequency'] = filtered_df['motif_frequency'].apply(lambda x: 1 if x >= bin_freq else 0)
    
    # Get all unique contigs that pass the filters (retain all contigs)
    all_contigs_info = filtered_df[['contig', 'contig_type', 'prefix', 'phylum', 'lineage', 'mge_type', 'depth', 'length']].drop_duplicates()
    all_mge_contigs = all_contigs_info[all_contigs_info['contig_type'] == 'MGE']['contig'].unique()
    all_host_contigs = all_contigs_info[all_contigs_info['contig_type'] == 'Host']['contig'].unique()
    
    ## Only consider motifs with frequency = 1 (present motifs)
    present_motifs = filtered_df[filtered_df['motif_frequency'] == 1]
    # print ("###", present_motifs)
    
    # Preserve contigs with no present motifs by adding placeholder rows
    contigs_with_motifs = set(present_motifs['contig'].unique())
    all_contigs = set(all_contigs_info['contig'].unique())
    contigs_without_motifs = all_contigs - contigs_with_motifs
    
    if len(contigs_without_motifs) > 0:
        # Create placeholder rows for contigs with no present motifs
        placeholder_rows = []
        for contig in contigs_without_motifs:
            contig_info = all_contigs_info[all_contigs_info['contig'] == contig].iloc[0]
            placeholder_rows.append({
                'motifString': None,
                'contig': contig,
                'motif_frequency': 0,
                'motif_loci_num': 0,
                'contig_type': contig_info['contig_type'],
                'mge_type': contig_info['mge_type'],
                'depth': contig_info['depth'],
                'length': contig_info['length'],
                'prefix': contig_info['prefix'],
                'phylum': contig_info['phylum'],
                'lineage': contig_info['lineage']
            })
        placeholder_df = pd.DataFrame(placeholder_rows)
        present_motifs = pd.concat([present_motifs, placeholder_df], ignore_index=True)
    
    mge_motifs = present_motifs[present_motifs['contig_type'] == 'MGE'][['motifString', 'contig', 'prefix', 
                                                                         'phylum', 'lineage','mge_type', 'depth', 'length','motif_loci_num']]
    host_motifs = present_motifs[present_motifs['contig_type'] == 'Host'][['motifString', 'contig', 'prefix', 
                                                                        'phylum', 'lineage','mge_type', 'depth', 'length','motif_loci_num']]

    # Use all MGE contigs that pass filters, not just those with present motifs
    filtered_mge_num = len(all_mge_contigs)
    
    jaccard_scores = []
    # Iterate over all MGE contigs, not just those with present motifs
    for mge_contig in all_mge_contigs:
        # Get motif set for this MGE (empty set if no present motifs)
        mge_motifs_set = set(mge_motifs[mge_motifs['contig'] == mge_contig]['motifString'])
        
        # Get host contig info
        if len(all_host_contigs) == 0:
            continue
        host_contig = all_host_contigs[0]
        host_motifs_set = set(host_motifs[host_motifs['contig'] == host_contig]['motifString'])
        intersection = mge_motifs_set.intersection(host_motifs_set)
        union = mge_motifs_set.union(host_motifs_set)
        # union = host_motifs_set ## only consider host motifs as MGE motifs might be error-prone
        ## exclude the motif that does not have motifstring in MGE
        # print (">>>",filtered_df)
        union_filtered = set()
        for motif in union:
            
            mge_motif_loci = filtered_df[(filtered_df['contig'] == mge_contig) & (filtered_df['motifString'] == motif)]['motif_loci_num']
            # print (mge_contig, "xxx", motif, mge_motif_loci)
            if len(mge_motif_loci) > 0 and mge_motif_loci.values[0] > 0:
                union_filtered.add(motif)

        if len(union) == 0:
            jaccard = 0.0
        else:
            jaccard = len(intersection) / len(union)

        if len(union_filtered) == 0:
            jaccard_filter = 0.0
        else:
            jaccard_filter = len(intersection) / len(union_filtered)
        relation = 'same_isolate'
        ## count number of motifs specifically in mge
        mge_specific = mge_motifs_set - host_motifs_set
        mge_specific_num = len(mge_specific)
        
        # Get contig metadata from all_contigs_info
        mge_info = all_contigs_info[all_contigs_info['contig'] == mge_contig].iloc[0]
        host_info = all_contigs_info[all_contigs_info['contig'] == host_contig].iloc[0]
        
        jaccard_scores.append({
            'mge_contig': mge_contig,
            'mge_type': mge_info['mge_type'],
            'host_contig': host_contig,
            'jaccard_similarity': jaccard,
            'jaccard_similarity_filtered': jaccard_filter,
            'share_num': len(intersection),
            'all_num': len(union),
            'all_num_filtered': len(union_filtered),
            'relation': relation,
            'mge_specific_num': mge_specific_num,
            'host_depth': host_info['depth'],
            'host_length': host_info['length'],
            'mge_depth': mge_info['depth'],
            'mge_length': mge_info['length'],
            'host_lineage': host_info['lineage'],
        })
    # print (jaccard_scores)
    ## count how many for each type of relation
    jaccard_scores = pd.DataFrame(jaccard_scores)

    return jaccard_scores, present_motifs, filtered_mge_num


def standardize_jaccard(mge_contig, host_contig, host_motifs, mge_motifs, bin_freq=0.3, 
                       motif_cache=None, profile_cache=None):
    """
    Optimized version with caching support.
    Pass motif_cache and profile_cache dicts to reuse across multiple calls.
    """
    if motif_cache is None:
        motif_cache = {}
    if profile_cache is None:
        profile_cache = {}
    
    host_motifs_set = set(host_motifs[host_motifs['contig'] == host_contig]['motifString'])
    mge_motifs_set = set(mge_motifs[mge_motifs['contig'] == mge_contig]['motifString'])
    union_motifs = mge_motifs_set.union(host_motifs_set)
    
    all_dir = "/home/shuaiw/borg/paper/isolation/batch2_results/"
    
    # Check cache for this contig pair
    pair_key = f"{mge_contig}_{host_contig}"
    if pair_key in profile_cache:
        re_profile_df = profile_cache[pair_key]
    else:
        # Build motif list efficiently using set for O(1) lookups
        motif_set_built = set()
        motif_list = []
        
        # Get prefixes once
        mge_prefix = mge_motifs[mge_motifs['contig'] == mge_contig]['prefix'].iloc[0]
        host_prefix = host_motifs[host_motifs['contig'] == host_contig]['prefix'].iloc[0]
        
        # Cache motif dataframes per contig
        if mge_contig not in motif_cache:
            ctg_obj = My_contig(mge_prefix, all_dir, mge_contig, data_type="isolation")
            motif_cache[mge_contig] = ctg_obj.read_motif()
        
        if host_contig not in motif_cache:
            ctg_obj = My_contig(host_prefix, all_dir, host_contig, data_type="isolation")
            motif_cache[host_contig] = ctg_obj.read_motif()
        
        # Build motif list from cached data
        for motif in union_motifs:
            if motif not in motif_set_built:
                # Check MGE motif data
                if motif in motif_cache[mge_contig]['motifString'].values:
                    centpos = motif_cache[mge_contig][motif_cache[mge_contig]['motifString'] == motif]['centerPos'].iloc[0]
                    motif_list.append([motif, centpos])
                    motif_set_built.add(motif)
                # Check host motif data
                elif motif in motif_cache[host_contig]['motifString'].values:
                    centpos = motif_cache[host_contig][motif_cache[host_contig]['motifString'] == motif]['centerPos'].iloc[0]
                    motif_list.append([motif, centpos])
                    motif_set_built.add(motif)
        
        prefix_list = [[mge_prefix, mge_contig], [host_prefix, host_contig]]
        
        # Call clean_profile and cache result
        re_profile_df = clean_profile(prefix_list, motif_list, all_dir, data_type="isolation")
        profile_cache[pair_key] = re_profile_df
    
    # Transform to binary based on bin_freq
    re_profile_df['binary_fraction'] = (re_profile_df['fraction'] >= bin_freq).astype(int)
    
    # Use boolean indexing more efficiently
    mge_mask = (re_profile_df['contig'] == mge_contig) & (re_profile_df['binary_fraction'] == 1)
    host_mask = (re_profile_df['contig'] == host_contig) & (re_profile_df['binary_fraction'] == 1)
    
    mge_motif_filter_set = set(re_profile_df[mge_mask]['motifString'].values)
    host_motif_filter_set = set(re_profile_df[host_mask]['motifString'].values)
    
    # Get jaccard similarity
    intersection = mge_motif_filter_set.intersection(host_motif_filter_set)
    union = mge_motif_filter_set.union(host_motif_filter_set)
    jaccard_similarity = len(intersection) / len(union) if len(union) > 0 else 0.0
    
    # Filter by motifs that actually exist in MGE
    mge_has_motif_mask = (re_profile_df['contig'] == mge_contig) & (re_profile_df['motif_loci_num'] > 0)
    mge_has_string_motif_set = set(re_profile_df[mge_has_motif_mask]['motifString'].values)
    
    # Remove motifs not in mge_has_string_motif_set
    mge_motif_filter_set = mge_motif_filter_set.intersection(mge_has_string_motif_set)
    host_motif_filter_set = host_motif_filter_set.intersection(mge_has_string_motif_set)
    
    # Get jaccard similarity after filtering
    intersection = mge_motif_filter_set.intersection(host_motif_filter_set)
    union = mge_motif_filter_set.union(host_motif_filter_set)
    jaccard_similarity_filtered = len(intersection) / len(union) if len(union) > 0 else 0.0
    # print ("###", intersection, union, jaccard_similarity_filtered, mge_has_motif_mask, mge_has_string_motif_set)
    # print (re_profile_df)
    return jaccard_similarity, jaccard_similarity_filtered


def process_mge_contig(mge_contig, mge_motifs, host_motifs, host_contig_list, 
                       clu_prefix_dict, random_ctg_num, bin_freq, max_per_relation=1):
    """Worker function to process one MGE contig against random host contigs."""
    # Each process has its own cache
    motif_cache = {}
    profile_cache = {}
    
    mge_motifs_set = set(mge_motifs[mge_motifs['contig'] == mge_contig]['motifString'])
    results = []
    test_count = 0
    
    # Track number of host contigs processed per relation
    relation_counts = {}
    
    for host_contig in np.random.choice(host_contig_list, size=min(random_ctg_num, 
                                                                   len(host_contig_list)), replace=False):
        
        # Get the prefix for the current MGE and host contigs
        mge_prefix = mge_motifs[mge_motifs['contig'] == mge_contig]['prefix'].iloc[0]
        host_prefix = host_motifs[host_motifs['contig'] == host_contig]['prefix'].iloc[0]
        
        if mge_prefix == host_prefix:
            relation = 'same_isolate'
        else:
            relation = esti_relation(
                mge_motifs[mge_motifs['contig'] == mge_contig]['lineage'].values[0],
                host_motifs[host_motifs['contig'] == host_contig]['lineage'].values[0]
            )
            if relation == "same_species":
                # check if they are in the same drep cluster
                if est_strain(clu_prefix_dict, mge_prefix, host_prefix):
                    test_count += 1
                    relation = 'same_strain'
        
        # Check if we've reached the max for this relation
        if relation_counts.get(relation, 0) >= max_per_relation:
            continue

        host_motifs_set = set(host_motifs[host_motifs['contig'] == host_contig]['motifString'])
        intersection = mge_motifs_set.intersection(host_motifs_set)
        union = mge_motifs_set.union(host_motifs_set)
        
        ## exclude the motif that does not have motifstring in MGE
        union_filtered = set()
        for motif in union:
            mge_motif_loci = mge_motifs[(mge_motifs['contig'] == mge_contig) & 
                                        (mge_motifs['motifString'] == motif)]['motif_loci_num']
            if len(mge_motif_loci) > 0 and mge_motif_loci.values[0] > 0:
                union_filtered.add(motif)

        jaccard = len(intersection) / len(union) if len(union) > 0 else 0.0
        jaccard_filter = len(intersection) / len(union_filtered) if len(union_filtered) > 0 else 0.0

        # Use cached version for speedup
        jaccard_similarity, jaccard_similarity_filtered = standardize_jaccard(
            mge_contig, host_contig, host_motifs, mge_motifs, bin_freq,
            motif_cache=motif_cache, profile_cache=profile_cache)
        
        # Increment relation counter
        relation_counts[relation] = relation_counts.get(relation, 0) + 1
        
        mge_type = mge_motifs[mge_motifs['contig'] == mge_contig]['mge_type'].values[0]
        results.append({
            'mge_contig': mge_contig,
            'mge_type': mge_type,
            'host_contig': host_contig,
            'jaccard_similarity': jaccard,
            'jaccard_similarity_filtered': jaccard_filter,
            'jaccard_similarity2': jaccard_similarity,
            'jaccard_similarity_filtered2': jaccard_similarity_filtered,
            'share_num': len(intersection),
            'all_num': len(union),
            'all_num_filtered': len(union_filtered),
            'relation': relation,
            'host_depth': host_motifs[host_motifs['contig'] == host_contig]['depth'].values[0],
            'host_length': host_motifs[host_motifs['contig'] == host_contig]['length'].values[0],
            'mge_depth': mge_motifs[mge_motifs['contig'] == mge_contig]['depth'].values[0],
            'mge_length': mge_motifs[mge_motifs['contig'] == mge_contig]['length'].values[0],
            'host_lineage': host_motifs[host_motifs['contig'] == host_contig]['lineage'].values[0],
        })
    
    return results, len(motif_cache), len(profile_cache)


def cross_sample_jaccard(present_motifs_all, clu_prefix_dict, random_ctg_num=100, bin_freq=0.3, n_workers=None):
    mge_motifs = present_motifs_all[present_motifs_all['contig_type'] == 'MGE'][['motifString', 'contig', 
                                                                                 'prefix', 'phylum', 'lineage','mge_type', 'depth', 'length', 'motif_loci_num']]
    host_motifs = present_motifs_all[present_motifs_all['contig_type'] == 'Host'][['motifString', 'contig', 
                                                                                   'prefix', 'phylum', 'lineage', 'mge_type', 'depth', 'length', 'motif_loci_num']]
    host_contig_list = host_motifs['contig'].unique()
    print(f"[✔] Found {len(mge_motifs)} MGE motifs and {len(host_motifs)} Host motifs with frequency = 1")
    print (f"[✔] Total unique host contigs: {len(host_contig_list)}")
    print (f"[✔] Total unique mge contigs: {len(mge_motifs['contig'].unique())}")
    ## count how many prefix in present_motifs_all
    print (f"[✔] Total unique prefixes: {len(present_motifs_all['prefix'].unique())}")

    # Use multiprocessing to parallelize MGE contig processing
    if n_workers is None:
        n_workers = min(os.cpu_count() or 4, len(mge_motifs['contig'].unique()))
    
    print(f"[✔] Using {n_workers} workers for parallel processing")
    
    mge_contig_list = mge_motifs['contig'].unique()
    total_pairs = len(mge_contig_list) * min(random_ctg_num, len(host_contig_list))
    processed = 0
    jaccard_scores = []
    total_motif_cache = 0
    total_profile_cache = 0
    
    # Use ProcessPoolExecutor for parallel processing
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        # Submit all MGE contigs for processing
        future_to_mge = {
            executor.submit(
                process_mge_contig, 
                mge_contig, mge_motifs, host_motifs, host_contig_list,
                clu_prefix_dict, random_ctg_num, bin_freq
            ): mge_contig for mge_contig in mge_contig_list
        }
        
        # Process results as they complete
        for future in as_completed(future_to_mge):
            mge_contig = future_to_mge[future]
            try:
                results, motif_cache_size, profile_cache_size = future.result()
                jaccard_scores.extend(results)
                processed += len(results)
                total_motif_cache += motif_cache_size
                total_profile_cache += profile_cache_size
                
                # Progress update
                if processed % 100 == 0 or processed == total_pairs:
                    print(f"[{processed}/{total_pairs}] Completed MGE contig {mge_contig}")
            except Exception as exc:
                print(f"[!] MGE contig {mge_contig} generated an exception: {exc}")
    
    print(f"\n[✔] Processed {processed} pairs total")
    print(f"[✔] Average motif cache size per worker: {total_motif_cache / len(mge_contig_list):.1f} contigs")
    print(f"[✔] Average profile cache size per worker: {total_profile_cache / len(mge_contig_list):.1f} pairs")
    
    jaccard_scores = pd.DataFrame(jaccard_scores)
    if len(jaccard_scores) > 0:
        relation_counts = jaccard_scores['relation'].value_counts()
        print(f"[✔] Relation counts:\n{relation_counts}")
    else:
        print("[!] No jaccard scores calculated (no MGE or host contigs found)")

    return jaccard_scores

def get_strain_prefix(drep_clu_dict):
    clu_prefix_dict = {}
    for clu_id in drep_clu_dict:
        clu_prefix_dict[clu_id] = set()
        for genome in drep_clu_dict[clu_id]:
            clu_prefix_dict[clu_id].add(genome.split("_")[0])
    return clu_prefix_dict
            
def est_strain(clu_prefix_dict, mge_prefix, host_prefix):
    for clu_id in clu_prefix_dict:
        if mge_prefix in clu_prefix_dict[clu_id] and host_prefix in clu_prefix_dict[clu_id]:
            return True
    return False

def esti_relation(lineage1, lineage2):
    """
    Estimate the relationship between two lineages based on taxonomic levels.
    
    Args:
        lineage1 (str): Taxonomic lineage of the first contig
        lineage2 (str): Taxonomic lineage of the second contig
    
    Returns:
        str: Estimated relationship (e.g., 'same_species', 'same_genus', etc.)
    """
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    ## either lineage is unknown, return different_lineage

    if re.search(r'Unclassified Bacteria', lineage1) or re.search(r'Unclassified Bacteria', lineage2):
        return 'different_lineage'

    taxa1 = lineage1.split(';')
    taxa2 = lineage2.split(';')
    # print ("###", lineage1, lineage2)
    
    for i in range(len(levels)-1, -1, -1):
        if taxa1[i] == taxa2[i] and taxa1[i] != 'unknown' and len(taxa1[i]) > 3:
            return f'same_{levels[i]}'
    
    return 'different_lineage'


def plot_motif_freq(all_motif_freq_df):
    ## plot clustered heatmap
    motif_pivot = all_motif_freq_df.pivot(index='prefix', columns='motif', values='frequency').fillna(0)
    
    # Create contig type annotation for rows
    # Get the contig type for each contig (prefix)
    contig_type_df = all_motif_freq_df[['prefix', 'contig_type']].drop_duplicates().set_index('prefix')
    
    # Create a color map for contig types
    contig_type_colors = {'Host': '#1f77b4', 'MGE': '#ff7f0e'}  # Blue for Host, Orange for MGE
    row_colors = contig_type_df['contig_type'].map(contig_type_colors)
    
    # Print data statistics to understand the range
    print(f"Data shape: {motif_pivot.shape}")
    print(f"Data range: {motif_pivot.min().min():.4f} to {motif_pivot.max().max():.4f}")
    print(f"Data mean: {motif_pivot.mean().mean():.4f}")
    print("Sample of data:")
    print(motif_pivot.head())
    print(f"Contig types: {contig_type_df['contig_type'].value_counts()}")
    
    # Clean the data to handle edge cases that cause NaN in clustering
    # Remove columns (motifs) that are all zeros
    motif_pivot = motif_pivot.loc[:, motif_pivot.sum(axis=0) > 0]
    
    # Remove rows (samples) that are all zeros
    motif_pivot = motif_pivot.loc[motif_pivot.sum(axis=1) > 0, :]
    
    # Align row colors with the cleaned data
    row_colors = row_colors.loc[motif_pivot.index]
    
    # Add small epsilon to avoid zero vectors (which cause issues with cosine distance)
    epsilon = 1e-10
    motif_pivot = motif_pivot + epsilon
    
    # Normalize by row (each sample sums to 1)
    # motif_pivot = motif_pivot.div(motif_pivot.sum(axis=1), axis=0)
    
    # Check for any remaining NaN values
    if motif_pivot.isna().any().any():
        print("Warning: Found NaN values after preprocessing, filling with zeros")
        motif_pivot = motif_pivot.fillna(0)
    
    print(f"After preprocessing - Data shape: {motif_pivot.shape}")
    print(f"After preprocessing - Data range: {motif_pivot.min().min():.6f} to {motif_pivot.max().max():.6f}")
    
    # Create a clustered heatmap without row standardization (since we already normalized)
    # Use euclidean distance instead of cosine to avoid numerical issues
    g = sns.clustermap(motif_pivot, cmap='viridis', metric='euclidean', method='average', 
                       standard_scale=None, figsize=(20, 20),  # No additional standardization
                       xticklabels=True, yticklabels=True,  # Ensure labels are shown
                       row_cluster=True, col_cluster=True,  # Keep clustering
                       dendrogram_ratio=(0.01, 0.01),  # Very small dendrograms
                       cbar_pos=(0.92, 0.3, 0.03, 0.4),  # Position colorbar outside: (x, y, width, height)
                       row_colors=row_colors,  # Add row annotation for contig types
                       colors_ratio=0.01)  # Minimize the width of the row color bar
    
    # # Hide the dendrograms by making them invisible
    # g.ax_row_dendrogram.set_visible(False)
    # g.ax_col_dendrogram.set_visible(False)
    
    # Rotate x-axis labels for better readability
    g.ax_heatmap.tick_params(axis='x', labelsize=12)
    g.ax_heatmap.tick_params(axis='y', rotation=45, labelsize=12)
    
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # Remove the colorbar/legend
    g.cax.set_visible(False)
    
    # Colorbar is now positioned outside via cbar_pos parameter
    # Add a label to the colorbar
    g.cax.set_ylabel('Motif Frequency', rotation=270, labelpad=20)
    
    # Create a custom legend for contig types
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=contig_type_colors['Host'], label='Host'),
                      Patch(facecolor=contig_type_colors['MGE'], label='MGE')]
    g.ax_heatmap.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))
    
    plt.title('Motif Frequency Heatmap')
    plt.savefig("../../tmp/results2/motif_freq_iso2.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def analyze_link(all_link_df):
    ## keep these with depth > 10
    # all_link_df = all_link_df[all_link_df['MGE_cov'] > 10]
    all_link_df = all_link_df[all_link_df['MGE_len'] > 10000]
    ## plot the distribution of final_score in all_link_df
    plt.figure(figsize=(10, 6))
    sns.histplot(all_link_df['final_score'], bins=50, kde=True, color='skyblue')
    plt.title('Distribution of Final Scores', fontsize=16, fontweight='bold')
    plt.xlabel('Final Score', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    plt.savefig("../../tmp/results2/final_score_distribution.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    ## output these with final_score == 0
    zero_score_df = all_link_df[all_link_df['final_score'] == 0]
    for  i, row in zero_score_df.iterrows():
        print (row["MGE"], row["host"], row["final_score"])
    
    ## three subplots, one is scatter plot of MGE_cov vs host_cov, one is MGE_gc vs host_gc, one is distribution of cos_sim
    print (len(all_link_df), "total links")
    fig, axs = plt.subplots(1, 3, figsize=(21, 7))
    # Scatter plot of MGE_cov vs host_cov
    axs[0].scatter(all_link_df['MGE_cov'], all_link_df['host_cov'], alpha=0.6, color='teal')
    axs[0].set_title('MGE Coverage vs Host Coverage', fontsize=16, fontweight='bold')
    axs[0].set_xlabel('MGE Coverage', fontsize=14)
    axs[0].set_ylabel('Host Coverage', fontsize=14)
    axs[0].grid(axis='both', alpha=0.3)
    # Set equal limits and add diagonal line for coverage plot
    cov_min = min(all_link_df['MGE_cov'].min(), all_link_df['host_cov'].min())
    cov_max = max(all_link_df['MGE_cov'].max(), all_link_df['host_cov'].max())
    axs[0].set_xlim(cov_min, cov_max)
    axs[0].set_ylim(cov_min, cov_max)
    axs[0].plot([cov_min, cov_max], [cov_min, cov_max], 'k--', alpha=0.7, linewidth=1.5, label='x = y')
    axs[0].legend()
    
    # Scatter plot of MGE_gc vs host_gc
    axs[1].scatter(all_link_df['MGE_gc'], all_link_df['host_gc'], alpha=0.6, color='coral')
    axs[1].set_title('MGE GC% vs Host GC%', fontsize=16, fontweight='bold')
    axs[1].set_xlabel('MGE GC%', fontsize=14)
    axs[1].set_ylabel('Host GC%', fontsize=14)
    axs[1].grid(axis='both', alpha=0.3)
    # Set equal limits and add diagonal line for GC plot
    gc_min = min(all_link_df['MGE_gc'].min(), all_link_df['host_gc'].min())
    gc_max = max(all_link_df['MGE_gc'].max(), all_link_df['host_gc'].max())
    axs[1].set_xlim(gc_min, gc_max)
    axs[1].set_ylim(gc_min, gc_max)
    axs[1].plot([gc_min, gc_max], [gc_min, gc_max], 'k--', alpha=0.7, linewidth=1.5, label='x = y')
    axs[1].legend()
    # Distribution of cos_sim
    sns.histplot(all_link_df['cos_sim'], bins=50, kde=True, color='orchid', ax=axs[2])
    axs[2].set_title('Distribution of Cosine Similarity', fontsize=16, fontweight='bold')
    axs[2].set_xlabel('Cosine Similarity', fontsize=14)
    axs[2].set_ylabel('Count', fontsize=14)
    axs[2].grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig("../../tmp/results2/link_features.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    for i, row in all_link_df.iterrows():
        if row["cos_sim"] < 0.5:
            print (row["MGE"], row["host"], row["MGE_gc"], row["host_gc"])

def plot_jaccard(jaccard_all, fig_dir):
    ## box plot hue to phylum
    plt.figure(figsize=(5, 5))
    sns.boxplot(data=jaccard_all, x='phylum', y='jaccard_similarity', hue='phylum', palette='Set3', legend=False)
    # plt.title('Jaccard Similarity between MGE and Host by Phylum', fontsize=16, fontweight='bold')
    plt.xlabel('Phylum', fontsize=14)
    plt.ylabel('Jaccard Similarity', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/jaccard_similarity_by_phylum.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    ## also a plot with two subplot, one for plasmid, and one for phage
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Subplot 1: Plasmid
    plasmid_data = jaccard_all[jaccard_all['mge_type'] == 'plasmid']
    if not plasmid_data.empty:
        sns.violinplot(data=plasmid_data, x='phylum', y='jaccard_similarity', 
                      ax=axes[0], palette='Set2')
        # Add mean markers
        for i, phylum in enumerate(plasmid_data['phylum'].unique()):
            mean_val = plasmid_data[plasmid_data['phylum'] == phylum]['jaccard_similarity'].mean()
            axes[0].scatter(i, mean_val, color='red', s=50, marker='D', zorder=10)
        axes[0].set_title('Plasmid')
        axes[0].set_xlabel('Phylum')
        axes[0].set_ylabel('Jaccard Similarity')
        axes[0].tick_params(axis='x', rotation=45)
        axes[0].grid(axis='y', alpha=0.3)
    
    # Subplot 2: Phage
    phage_data = jaccard_all[jaccard_all['mge_type'] == 'virus']
    if not phage_data.empty:
        sns.violinplot(data=phage_data, x='phylum', y='jaccard_similarity', 
                      ax=axes[1], palette='Set2')
        # Add mean markers
        for i, phylum in enumerate(phage_data['phylum'].unique()):
            mean_val = phage_data[phage_data['phylum'] == phylum]['jaccard_similarity'].mean()
            axes[1].scatter(i, mean_val, color='red', s=50, marker='D', zorder=10)
        axes[1].set_title('Virus')
        axes[1].set_xlabel('Phylum')
        axes[1].set_ylabel('Jaccard Similarity')
        axes[1].tick_params(axis='x', rotation=45)
        axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/jaccard_similarity_by_phylum_mge_type.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def gradient_plot(jaccard_all, fig_dir):
    ## plot the proportion of MGE-host pairs with at least 0, 1, 2, ... shared motifs  using bar plot
    max_share = jaccard_all['share_num'].max()
    share_counts = []
    for i in range(0, int(max_share) + 1):
        count = len(jaccard_all[jaccard_all['share_num'] >= i])
        share_counts.append({'shared_motifs': i, 'count': count, 'proportion': count / len(jaccard_all)})
    share_counts_df = pd.DataFrame(share_counts)
    plt.figure(figsize=(5, 5))
    sns.barplot(data=share_counts_df, x='shared_motifs', y='proportion', color='steelblue')
    # plt.title('Number of MGE-Host Pairs by Shared Motifs', fontsize=16, fontweight='bold')
    plt.xlabel('Number of Shared Motifs', fontsize=14)
    plt.ylabel('Proportion of MGE-Host Pairs', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/mge_host_shared_motifs.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    ## also a plot with two subplot, one for plasmid, and one for phage
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Subplot 1: Plasmid
    plasmid_data = jaccard_all[jaccard_all['mge_type'] == 'plasmid']
    if not plasmid_data.empty:
        plasmid_share_counts = []
        for i in range(0, int(max_share) + 1):
            count = len(plasmid_data[plasmid_data['share_num'] >= i])
            plasmid_share_counts.append({'shared_motifs': i, 'count': count, 'proportion': count / len(plasmid_data)})
        plasmid_share_counts_df = pd.DataFrame(plasmid_share_counts)
        sns.barplot(data=plasmid_share_counts_df, x='shared_motifs', y='proportion', color='lightcoral', ax=axes[0])
        axes[0].set_title('Plasmid')
        axes[0].set_xlabel('Number of Shared Motifs')
        axes[0].set_ylabel('Proportion of MGE-Host Pairs')
        axes[0].grid(axis='y', alpha=0.3)
    # Subplot 2: Phage
    phage_data = jaccard_all[jaccard_all['mge_type'] == 'virus']
    if not phage_data.empty:
        phage_share_counts = []
        for i in range(0, int(max_share) + 1):
            count = len(phage_data[phage_data['share_num'] >= i])
            phage_share_counts.append({'shared_motifs': i, 'count': count, 'proportion': count / len(phage_data)})
        phage_share_counts_df = pd.DataFrame(phage_share_counts)
        sns.barplot(data=phage_share_counts_df, x='shared_motifs', y='proportion', color='lightseagreen', ax=axes[1])
        axes[1].set_title('Virus')
        axes[1].set_xlabel('Number of Shared Motifs')
        axes[1].set_ylabel('Proportion of MGE-Host Pairs')
        axes[1].grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/mge_host_shared_motifs_by_mge_type.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def cross_taxa_plot(jaccard_all, fig_dir):
    ## box plot hue to relation with ordered x-axis
    # Define the desired order from most closely related to least related

    # Filter out rows with invalid relation values (NaN, None, or non-string)
    jaccard_all = jaccard_all[jaccard_all['relation'].notna()]
    jaccard_all = jaccard_all[jaccard_all['relation'].isin(relation_order)]
    
    # Create a figure with 2x2 subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    print ("Total pairs:", len(jaccard_all))
    # Top left subplot: All data
    sns.boxplot(data=jaccard_all, x='relation', y='jaccard_similarity', 
                hue='relation', palette='Set2', legend=False, order=relation_order, ax=ax1)
    ax1.set_xlabel('Taxonomic Relation', fontsize=14)
    ax1.set_ylabel('Jaccard Similarity', fontsize=14)
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_title('All Data', fontsize=14)
    
    # Top right subplot: Filtered data (share_num > 1)
    filtered_data = jaccard_all[jaccard_all['all_num'] > 2]
    print (len(filtered_data), "pairs with share union num > 2")
    sns.boxplot(data=filtered_data, x='relation', y='jaccard_similarity', 
                hue='relation', palette='Set2', legend=False, order=relation_order, ax=ax2)
    ax2.set_xlabel('Taxonomic Relation', fontsize=14)
    ax2.set_ylabel('Jaccard Similarity', fontsize=14)
    ax2.tick_params(axis='x', rotation=45)
    ax2.grid(axis='y', alpha=0.3)
    ax2.set_title('Union motif Num > 2', fontsize=14)
    
    # Bottom left subplot: All data violin plot
    sns.violinplot(data=jaccard_all, x='relation', y='jaccard_similarity', 
                  hue='relation', ax=ax3, palette='Set3', order=relation_order, legend=False)
    # Add mean markers
    for i, relation in enumerate(relation_order):
        relation_data = jaccard_all[jaccard_all['relation'] == relation]
        if not relation_data.empty:
            mean_val = relation_data['jaccard_similarity'].mean()
            ax3.scatter(i, mean_val, color='red', s=50, marker='D', zorder=10)
    ax3.set_title('All Data (Violin)', fontsize=14)
    ax3.set_xlabel('Taxonomic Relation', fontsize=14)
    ax3.set_ylabel('Jaccard Similarity', fontsize=14)
    ax3.tick_params(axis='x', rotation=45)
    ax3.grid(axis='y', alpha=0.3)
    
    # Bottom right subplot: Filtered data violin plot
    sns.violinplot(data=filtered_data, x='relation', y='jaccard_similarity', 
                  hue='relation', ax=ax4, palette='Set3', order=relation_order, legend=False)
    # Add mean markers
    for i, relation in enumerate(relation_order):
        relation_data = filtered_data[filtered_data['relation'] == relation]
        if not relation_data.empty:
            mean_val = relation_data['jaccard_similarity'].mean()
            ax4.scatter(i, mean_val, color='red', s=50, marker='D', zorder=10)
    ax4.set_title('Share motif Num > 1 (Violin)', fontsize=14)
    ax4.set_xlabel('Taxonomic Relation', fontsize=14)
    ax4.set_ylabel('Jaccard Similarity', fontsize=14)
    ax4.tick_params(axis='x', rotation=45)
    ax4.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/jaccard_similarity_by_taxa_relation.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    ## also plot a violin plot with two subplot, one for plasmid, and one for phage
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # Subplot 1: Plasmid
    plasmid_data = jaccard_all[jaccard_all['mge_type'] == 'plasmid']
    if not plasmid_data.empty:      
        sns.violinplot(data=plasmid_data, x='relation', y='jaccard_similarity', 
                      hue='relation', ax=axes[0], palette='Set3', order=relation_order, legend=False)
        # Add mean markers
        for i, relation in enumerate(relation_order):
            mean_val = plasmid_data[plasmid_data['relation'] == relation]['jaccard_similarity'].mean()
            axes[0].scatter(i, mean_val, color='red', s=50, marker='D', zorder=10)
        axes[0].set_title('Plasmid')
        axes[0].set_xlabel('Taxonomic Relation')
        axes[0].set_ylabel('Jaccard Similarity')
        axes[0].tick_params(axis='x', rotation=45)
        axes[0].grid(axis='y', alpha=0.3)
    # Subplot 2: Phage
    phage_data = jaccard_all[jaccard_all['mge_type'] == 'virus']
    if not phage_data.empty:    
        sns.violinplot(data=phage_data, x='relation', y='jaccard_similarity', 
                      hue='relation', ax=axes[1], palette='Set3', order=relation_order, legend=False)
        # Add mean markers
        for i, relation in enumerate(relation_order):
            mean_val = phage_data[phage_data['relation'] == relation]['jaccard_similarity'].mean()
            axes[1].scatter(i, mean_val, color='red', s=50, marker='D', zorder=10)
        axes[1].set_title('Virus')
        axes[1].set_xlabel('Taxonomic Relation')
        axes[1].set_ylabel('Jaccard Similarity')
        axes[1].tick_params(axis='x', rotation=45)
        axes[1].grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/jaccard_similarity_by_taxa_relation_violin.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def count_jaccard(same_sample_df, jaccard_all_sample):
    ## print all rows with jaccard similarity < 0.5 for check
    low_jaccard = same_sample_df[same_sample_df['jaccard_similarity'] < 0.5]
    print("Low Jaccard Similarity Samples:")
    for index, row in low_jaccard.iterrows():
        print(f"Prefix: {row['prefix']}, MGE Contig: {row['mge_contig']}, Host Contig: {row['host_contig']}, Jaccard Similarity: {row['jaccard_similarity']:.4f}")
    
    ## count proportion of MGE-host pairs with jaccard similarity =1
    jaccard_1 = same_sample_df[same_sample_df['jaccard_similarity'] == 1]
    print(f"Proportion of MGE-host pairs with Jaccard similarity = 1: {len(jaccard_1) / len(same_sample_df):.4f}")
    ## compare cross-sample and same-sample jaccard similarity using boxplot, jaccard_all_sample have both same_sample and cross_sample
    combined_df = pd.concat([same_sample_df.assign(sample_type='same_isolate'),
                             jaccard_all_sample[jaccard_all_sample['relation'] != 'same_isolate'].assign(sample_type='cross_sample')])
    plot_compare(combined_df, fig_dir)
    ## print the rows with mge_specific_num > 0
    mge_specific = same_sample_df[same_sample_df['mge_specific_num'] > 0]
    print("MGE-specific motifs in same-sample pairs:")
    for index, row in mge_specific.iterrows():
        print(f"Prefix: {row['prefix']}, MGE Contig: {row['mge_contig']}, Host Contig: {row['host_contig']}, MGE-specific Motifs: {row['mge_specific_num']}")
    ## count the number of MGE-host pairs with for each relation in jaccard_all_sample
    relation_counts = jaccard_all_sample['relation'].value_counts()
    print("Number of MGE-host pairs by taxonomic relation:")
    for relation, count in relation_counts.items():
        print(f"{relation}: {count}")

def plot_compare(combined_df, fig_dir):
    plt.figure(figsize=(5, 5))
    sns.boxplot(data=combined_df, x='sample_type', y='jaccard_similarity', hue='sample_type', palette='Set1', legend=False)
    plt.xlabel('Sample Type', fontsize=14)
    plt.ylabel('Jaccard Similarity', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    ## add statistical annotation
    from scipy.stats import mannwhitneyu
    same_sample_values = combined_df[combined_df['sample_type'] == 'same_isolate']['jaccard_similarity']
    cross_sample_values = combined_df[combined_df['sample_type'] == 'cross_sample']['jaccard_similarity']
    stat, p_value = mannwhitneyu(same_sample_values, cross_sample_values, alternative='two-sided')
    ## add p-value to the plot
    plt.text(0.5, max(combined_df['jaccard_similarity']) * 0.95, f'p-value = {p_value:.4e}', ha='center', fontsize=12)
    plt.savefig(f"{fig_dir}/jaccard_similarity_same_vs_cross_sample.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def plot_jaccard_distribution(same_sample_df, fig_dir):
    ## plot the distribution of jaccard similarity using histogram
    plt.figure(figsize=(5, 5))
    sns.histplot(same_sample_df['jaccard_similarity'], bins=20, kde=True, color='lightgreen')
    # plt.title('Distribution of Jaccard Similarity (Same Sample)', fontsize=16, fontweight='bold')
    plt.xlabel('Jaccard Similarity', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/jaccard_similarity_distribution_same_sample.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    ## a plot with two subplots, one for plasmid, and one for virus
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    # Subplot 1: Plasmid
    plasmid_data = same_sample_df[same_sample_df['mge_type'] == 'plasmid']
    if not plasmid_data.empty:
        sns.histplot(plasmid_data['jaccard_similarity'], bins=20, kde=True, color='salmon', ax=axes[0])
        axes[0].set_title('Plasmid')
        axes[0].set_xlabel('Jaccard Similarity')
        axes[0].set_ylabel('Count')
        axes[0].grid(axis='y', alpha=0.3)
    # Subplot 2: Virus
    virus_data = same_sample_df[same_sample_df['mge_type'] == 'virus']
    if not virus_data.empty:
        sns.histplot(virus_data['jaccard_similarity'], bins=20, kde=True, color='skyblue', ax=axes[1])
        axes[1].set_title('Virus')
        axes[1].set_xlabel('Jaccard Similarity')
        axes[1].set_ylabel('Count')
        axes[1].grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/jaccard_similarity_distribution_by_mge_type.pdf", dpi=300, bbox_inches="tight")
    plt.close()

def plot_MTase(df_all_data, fig_dir):
    ##scatter plot and box plot for MTase num vs motif num
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Scatter plot with regression line
    sns.scatterplot(data=df_all_data, x='RM_num', y='motif_num', 
                   palette='Set2', alpha=0.6, s=40, ax=ax1)
    # Add regression line
    sns.regplot(data=df_all_data, x='RM_num', y='motif_num', scatter=False, 
                color='red', ax=ax1)
    ax1.set_title('MTase Number vs Motif Number (Scatter)')
    ax1.set_xlabel('MTase Number')
    ax1.set_ylabel('Motif Number')
    ax1.legend(title='x', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Box plot - each x shows a number of MTase, y shows number of motifs
    # Sort the RM_num values and create ordered categories
    rm_order = sorted(df_all_data['RM_num'].unique())
    df_all_data['RM_num_str'] = df_all_data['RM_num'].astype(str)
    sns.boxplot(data=df_all_data, x='RM_num_str', y='motif_num',
                order=[str(x) for x in rm_order], ax=ax2)
    ax2.set_title('Motif Number Distribution by MTase Count')
    ax2.set_xlabel('MTase Number')
    ax2.set_ylabel('Motif Number')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/MTase_vs_motif_num.png", dpi=300, bbox_inches='tight')
    plt.close()

def plot_motif_num(df_all_data, fig_dir):
    ## box plot, x is phylum, y is motif num
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=df_all_data, x='phylum', y='motif_num', palette='Set3')
    plt.xlabel('Phylum', fontsize=14)
    plt.ylabel('Motif Number', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.savefig(f"{fig_dir}/motif_num_by_phylum.pdf", dpi=300, bbox_inches="tight")
    plt.close()


def main_profile(all_dir, fig_dir):
    """
    samples: bacteria, pure, high average depth >= 10x, has MGEs, has motifs
    """
    i = 0
    pure_num, mge_num, has_motif_num = 0, 0, 0
    jaccard_all = pd.DataFrame()
    present_motifs_all = pd.DataFrame()
    archea_list = ["SRR27457941", "SRR31014709"]
    drep_clu_dict = read_drep_cluster(drep_clu_file, {})
    clu_prefix_dict = get_strain_prefix(drep_clu_dict)
    profile_data = []
    potential_megaP_all = []
    # print (clu_prefix_dict)
    # print (est_strain(clu_prefix_dict, 'SRR32364911', 'SRR32364914'))
    
    for prefix in os.listdir(all_dir):
        if prefix in archea_list:
            continue
        i += 1
        sample_obj = Isolation_sample(prefix, all_dir)


        sample_obj.get_phylum()
        mge_dict = sample_obj.read_MGE()
        depth_dict, length_dict = sample_obj.read_depth()
        pure_flag = sample_obj.check_pure2()
        average_dp = sample_obj.get_average_depth()
        unique_motif_num, unique_motifs = sample_obj.get_unique_motifs()
        potential_megaP = sample_obj.search_megaP() 
        potential_megaP_all += potential_megaP 
        print (potential_megaP)
        
        # sample_obj.explore_specific_motifs()
        sample_obj.load_isolation_RM()
        MT_num = len(sample_obj.isolation_RM_dict)

        if pure_flag == "pure":
            pure_num += 1
        else:
            continue
        if average_dp < 10:
            print(f"Skipping {prefix} as low average depth: {average_dp}")
            continue
        # if len(mge_dict) < 1:
        #     continue
        # else:
        #     mge_num += 1

        if unique_motifs == 0:
            print(f"Skipping {prefix} as no motifs.")
            continue
        else:
            has_motif_num += 1

        if not os.path.exists(sample_obj.profile):
            continue
        line_count = sum(1 for line in open(sample_obj.profile)) - 1  # Subtract 1 for header
        if line_count < 2:
            print(f"Skipping {prefix} as insufficient motif data.")
            continue
        
        profile_data.append([prefix, sample_obj.phylum, unique_motif_num, MT_num])
    profile_df = pd.DataFrame(profile_data, columns=['sample', 'phylum', 'motif_num', 'RM_num'])
    profile_df.to_csv(f"{fig_dir}/isolation_samples_motif_profile.csv", index=False)
    

    print (len(potential_megaP_all), "potential megaP contigs found:")
    for contig_info in potential_megaP_all:
        print (contig_info[0], contig_info[1], contig_info[2], contig_info[3])

    # profile_df = pd.read_csv(f"{fig_dir}/isolation_samples_motif_profile.csv")
    # plot_MTase(profile_df, fig_dir)
    # plot_motif_num(profile_df, fig_dir)
   
def main(all_dir, fig_dir, drep_clu_file):
    """
    samples: bacteria, pure, high average depth >= 10x, has MGEs, has motifs
    """
    i = 0
    pure_num, mge_num, has_motif_num, final_num = 0, 0, 0, 0
    jaccard_all = pd.DataFrame()
    present_motifs_all = pd.DataFrame()
    all_contig_info = pd.DataFrame()
    archea_list = ["SRR27457941", "SRR31014709"]
    drep_clu_dict = read_drep_cluster(drep_clu_file, {})
    clu_prefix_dict = get_strain_prefix(drep_clu_dict)
    profile_data = []
    # print (clu_prefix_dict)
    # print (est_strain(clu_prefix_dict, 'SRR32364911', 'SRR32364914'))
    all_mge_dict = {}
    all_filtered_mge_num = 0

    for prefix in os.listdir(all_dir):
        if prefix in archea_list:
            continue
        # if prefix != "SRR29842206":
        #     continue
        i += 1
        sample_obj = Isolation_sample(prefix, all_dir)

        sample_obj.get_phylum()
        mge_dict = sample_obj.read_MGE()
        depth_dict, length_dict = sample_obj.read_depth()
        pure_flag = sample_obj.check_pure2()
        average_dp = sample_obj.get_average_depth()
        unique_motif_num, unique_motifs = sample_obj.get_unique_motifs()
        sample_obj.explore_specific_motifs()
        genome_list, represent_contig_list = sample_obj.get_iso_good_ctgs(min_depth=10, min_len=500000)

        if average_dp < 10:
            print(f"Skipping {prefix} as low average depth: {average_dp}")
            continue

        if pure_flag == "pure":
            pure_num += 1
        else:
            continue

        # print (f"{prefix} is {pure_flag}")
        if len(mge_dict) < 1:
            continue
        else:
            mge_num += 1

        if len(unique_motifs) == 0:
            print(f"Skipping {prefix} as no motifs.")
            continue
        else:
            has_motif_num += 1

        if len(represent_contig_list) < 1:
            print(f"Skipping {prefix} as no good contigs.")
            continue

        # Enhanced motif sharing analysis with additional annotations
        jaccard_scores, present_motifs, annotations_df, filtered_mge_num = quantify_sharing(sample_obj, mge_dict, depth_dict, 
                                                          length_dict, unique_motifs, 
                                                          represent_contig_list)
        # print (annotations_df)
        ## check mge number in filtered_df, skip if zero
        filter_mge_num = len(annotations_df[annotations_df['contig_type'] == "MGE"])
        if filter_mge_num == 0:
            print(f"Skipping {prefix} as no MGEs in filtered dataframe.")
            continue
        all_filtered_mge_num += filtered_mge_num
        jaccard_scores['prefix'] = prefix
        jaccard_scores['phylum'] = sample_obj.phylum
        all_contig_info = pd.concat([all_contig_info, annotations_df], axis=0, ignore_index=True)
        all_mge_dict = {**all_mge_dict, **mge_dict}
        # Use robust concatenation that handles empty DataFrames
        if jaccard_all.empty:
            jaccard_all = jaccard_scores.copy()
        else:
            jaccard_all = pd.concat([jaccard_all, jaccard_scores], axis=0, ignore_index=True)
            
        if present_motifs_all.empty:
            present_motifs_all = present_motifs.copy()
        else:
            present_motifs_all = pd.concat([present_motifs_all, present_motifs], axis=0, ignore_index=True)
        
        final_num += 1

        # if i > 50:
        #     break  # For testing, limit to first 50 samples

    # count_mge(all_mge_dict, all_contig_info)
    
    same_sample_df = jaccard_all[jaccard_all['relation'] == 'same_isolate']
    print (same_sample_df)
    print (len(same_sample_df) , "total same-isolate MGE-host pairs out of ", len(jaccard_all))
    
    jaccard_all_sample = cross_sample_jaccard(present_motifs_all, clu_prefix_dict, random_ctg_num=1000, bin_freq=0.3, n_workers=32)
    jaccard_all_sample.to_csv(f"{fig_dir}/jaccard_all_samples.csv", index=False)
    same_sample_df.to_csv(f"{fig_dir}/jaccard_same_sample.csv", index=False)


    # jaccard_all_sample = pd.read_csv(f"{fig_dir}/jaccard_all_samples.csv")
    # same_sample_df = pd.read_csv(f"{fig_dir}/jaccard_same_sample.csv")
    # # analyze_link(all_link_df)
    # plot_jaccard(same_sample_df, fig_dir)
    # gradient_plot(same_sample_df, fig_dir)
    # plot_jaccard_distribution(same_sample_df, fig_dir)

    # cross_taxa_plot(jaccard_all_sample, fig_dir)
    # # count_jaccard(same_sample_df, jaccard_all_sample)
    # print (f"Total {i} samples, {pure_num} pure samples, {mge_num} samples with MGEs, {has_motif_num} samples with motifs, {final_num} samples analyzed.")
    # print (f"Total filtered MGEs: {all_filtered_mge_num}")
    # print ("Analysis complete.")


def get_new_host(plasmid_list):
    plasmid_host_dict = {}
    df = pd.read_csv(plasmid_list)
    for index, row in df.iterrows():
        plasmid = row['seq_name']
        host_str = row['host']

        plasmid_host_dict[plasmid] = ['', '', host_str.split(";")]
    return plasmid_host_dict


def main_96plex(fig_dir, bin_freq=0.3):
    data_dir ="/home/shuaiw/borg/paper/linkage/pure2/m64004_210929_143746.p100/"
    plasmid_list_file = "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list"
    my_96plex_data = []
    plasmid_host_dict = get_new_host(plasmid_list_file)
    motif_profile = os.path.join(data_dir, "motif_profile.csv")
    profile_df = pd.read_csv(motif_profile)
    
    # Melt the profile_df from wide to long format
    profile_df_melted = profile_df.melt(id_vars=['motif_identifier'], 
                                        var_name='contig', 
                                        value_name='freq')
    
    ## binary profile_df with freq >= 0.3 as 1, else 0
    profile_df_melted = profile_df_melted[profile_df_melted['freq'] >= bin_freq]

    for my_MGE in plasmid_host_dict:
        host = plasmid_host_dict[my_MGE][2][0]
        print (my_MGE, host)
        my_MGE_motif_file = os.path.join(data_dir, f"motifs/{my_MGE}.motifs.csv")
        my_host_motif_file = os.path.join(data_dir, f"motifs/{host}.motifs.csv")
        if not os.path.exists(my_MGE_motif_file) or not os.path.exists(my_host_motif_file):
            continue

        my_host_motif_df = pd.read_csv(my_host_motif_file)
        host_motif_num, host_motifs, unique_motifs_identifier = get_unique_motifs(my_host_motif_df)

        ## get motif set for MGE and host, based on profile_df_melted,
        # # only consider motifs in host_motifs
        MGE_motif_set = set(profile_df_melted[
            (profile_df_melted['contig'] == my_MGE) &
            (profile_df_melted['motif_identifier'].isin(unique_motifs_identifier))
        ]['motif_identifier'])
        host_motif_set = set(profile_df_melted[
            (profile_df_melted['contig'] == host) &
            (profile_df_melted['motif_identifier'].isin(unique_motifs_identifier))
        ]['motif_identifier'])

        intersection = MGE_motif_set.intersection(host_motif_set)
        union = MGE_motif_set.union(host_motif_set)
        if len(union) == 0:
            jaccard = 0.0
        else:
            jaccard = len(intersection) / len(union)
        my_96plex_data.append([my_MGE, host, len(intersection), len(union), jaccard])
    my_96plex_df = pd.DataFrame(my_96plex_data, columns=['MGE', 'host', 'shared_motif_num', 'total_motif_num', 'jaccard_similarity'])
    print (my_96plex_df)
    plot_96plex(my_96plex_df, fig_dir)
    

def plot_96plex(my_96plex_df, fig_dir):
    ## plot the distribution of jaccard similarity using histogram
    # Create subplot with histogram and barplot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot histogram of jaccard similarity
    sns.histplot(my_96plex_df['jaccard_similarity'], bins=20, kde=True, color='lightcoral', ax=ax1)
    ax1.set_xlabel('Jaccard Similarity', fontsize=14)
    ax1.set_ylabel('Count', fontsize=14)
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_title('Distribution of Jaccard Similarity', fontsize=14)
    
    # Plot barplot of shared motif number cutoff
    max_shared = my_96plex_df['shared_motif_num'].max()
    
    # Check if max_shared is NaN or if dataframe is empty
    if pd.isna(max_shared) or len(my_96plex_df) == 0:
        print("Warning: No valid shared motif data found, skipping barplot")
        ax2.text(0.5, 0.5, 'No Data Available', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_xlabel('Number of Shared Motifs', fontsize=14)
        ax2.set_ylabel('Proportion of MGE-Host Pairs', fontsize=14)
    else:
        share_counts = []
        for i in range(0, int(max_shared) + 1):
            count = len(my_96plex_df[my_96plex_df['shared_motif_num'] >= i])
            share_counts.append({'shared_motifs': i, 'count': count, 'proportion': count / len(my_96plex_df)})
        share_counts_df = pd.DataFrame(share_counts)
        
        sns.barplot(data=share_counts_df, x='shared_motifs', y='proportion', color='steelblue', ax=ax2)
        ax2.set_xlabel('Number of Shared Motifs', fontsize=14)
        ax2.set_ylabel('Proportion of MGE-Host Pairs', fontsize=14)
        ax2.grid(axis='y', alpha=0.3)
        ax2.set_title('Shared Motifs Distribution', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/jaccard_similarity_distribution_96plex.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    my_96plex_df.to_csv(f"{fig_dir}/jaccard_similarity_96plex.csv", index=False)

def analyze_jaccard(fig_dir):
    jaccard_all_sample = pd.read_csv(f"{fig_dir}/jaccard_all_samples.csv")
    same_sample_df = pd.read_csv(f"{fig_dir}/jaccard_same_sample.csv")  


    # same_sample_df = same_sample_df[same_sample_df['mge_length'] >= 50000]


    ## calculate the proportion of rows with jaccard_similarity==1 in same_sample_df
    perfect_match = same_sample_df[same_sample_df['jaccard_similarity_filtered'] == 1]
    proportion_perfect = len(perfect_match) / len(same_sample_df)
    print(f"\nProportion of MGE-host pairs with perfect Jaccard similarity (=1): {proportion_perfect:.4f} ({len(perfect_match)}/{len(same_sample_df)})")
    ## get the perfect_match proportion for each relation  in jaccard_all_sample
    relation_groups = jaccard_all_sample.groupby('relation')
    print("\nProportion of perfect Jaccard similarity (=1) by taxonomic relation:")
    data = []
    for relation, group in relation_groups:
        perfect_count = len(group[group['jaccard_similarity_filtered'] == 1])
        total_count = len(group)
        proportion = perfect_count / total_count if total_count > 0 else 0
        print(f"{relation}: {proportion:.4f} ({perfect_count}/{total_count})")
        data.append([relation, perfect_count, total_count, proportion])
    relation_df = pd.DataFrame(data, columns=['relation', 'perfect_count', 'total_count', 'proportion'])
    ## plot bar plot for relation_df
    plt.figure(figsize=(8, 6))
    sns.barplot(data=relation_df, x='relation', y='proportion', palette='Set2', order=relation_order)
    plt.xlabel('Taxonomic Relation', fontsize=14)
    plt.ylabel('Proportion of Perfect Jaccard Similarity (=1)', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.title('Proportion of Perfect Jaccard Similarity by Taxonomic Relation', fontsize=16, fontweight='bold')
    plt.savefig(f"{fig_dir}/proportion_perfect_jaccard_by_relation.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    ## get the perfect_match proportion for each phylum in same_sample_df
    phylum_groups = same_sample_df.groupby('phylum')
    print("\nProportion of perfect Jaccard similarity (=1) by phylum:")
    data = []
    for phylum, group in phylum_groups:
        perfect_count = len(group[group['jaccard_similarity_filtered'] == 1])
        total_count = len(group)
        proportion = perfect_count / total_count if total_count > 0 else 0
        print(f"{phylum}: {proportion:.4f} ({perfect_count}/{total_count})")
        data.append([phylum, perfect_count, total_count, proportion])
    phylum_df = pd.DataFrame(data, columns=['phylum', 'perfect_count', 'total_count', 'proportion'])
    ## plot bar plot for phylum_df
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(data=phylum_df, x='phylum', y='proportion', palette='Set3')
    
    # Add total_count labels on top of each bar
    for i, (idx, row) in enumerate(phylum_df.iterrows()):
        ax.text(i, row['proportion'] + 0.02, f"n={int(row['total_count'])}", 
                ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.xlabel('Phylum', fontsize=14)
    plt.ylabel('Proportion of Perfect Jaccard Similarity (=1)', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    plt.title('Proportion of Perfect Jaccard Similarity by Phylum', fontsize=16, fontweight='bold')
    plt.savefig(f"{fig_dir}/proportion_perfect_jaccard_by_phylum.pdf", dpi=300, bbox_inches="tight")
    plt.close()

    # ## see the rows that jaccard_similarity < 0.6 in same_sample_df
    # low_jaccard = same_sample_df[same_sample_df['jaccard_similarity'] < 1]
    # for index, row in low_jaccard.iterrows():
    #     if row['phylum'] != 'Bacillota':
    #         continue
    #     print(f"\n--- Row {index} ---")
    #     for col in row.index:
    #         print(f"{col}: {row[col]}")


def count_mge(all_mge_dict, all_contig_info):
    ## count the mge types in all_mge_dict
    mge_type_counts = {}
    for mge_contig, mge_info in all_mge_dict.items():
        mge_type = mge_info[0]
        if mge_type not in mge_type_counts:
            mge_type_counts[mge_type] = 0
        mge_type_counts[mge_type] += 1
    print("MGE type counts across all samples:", len(all_mge_dict))
    for mge_type, count in mge_type_counts.items():
        print(f"{mge_type}: {count}")
        print(f"Proportion of {mge_type}: {count/len(all_mge_dict):.4f}")

    ## count number of host contigs
    host_contigs = all_contig_info[all_contig_info['contig_type'] == "Host"]
    unique_host_contigs = host_contigs['contig'].nunique()
    print(f"\nTotal unique host contigs: {unique_host_contigs}")
    ## count the mge types in all_contig_info
    all_contig_info_mge = all_contig_info[all_contig_info['contig_type'] == "MGE"]
    print(f"\nTotal MGE contigs in all_contig_info: {len(all_contig_info_mge)}")
    print(f"Total unique MGE contigs in all_contig_info: {all_contig_info_mge['contig'].nunique()}")
    contig_mge_type_counts = all_contig_info_mge['mge_type'].value_counts()
    print("\nMGE type counts in filtered contig info:")
    for mge_type, count in contig_mge_type_counts.items():
        print(f"{mge_type}: {count}")
        print(f"Proportion of {mge_type}: {count/len(all_contig_info_mge):.4f}")
    

if __name__ == "__main__":
    # meta_file = "/home/shuaiw/Methy/assembly_pipe/prefix_table.tab"
    # sample_env_dict = read_metadata(meta_file)
    # all_dir = "/home/shuaiw/borg/paper/isolation/bacteria/"
    all_dir = "/home/shuaiw/borg/paper/isolation//batch2_results/"
    fig_dir = "../../tmp/figures/motif_sharing/"
    profile_fig_dir = "../../tmp/figures/iso_profile/"

    ANI=99
    drep_clu_file = f"/home/shuaiw/borg/paper/specificity/iso_{ANI}_out/data_tables/Cdb.csv"

    # main(all_dir, fig_dir, drep_clu_file)
    main_96plex(fig_dir)
    
    # main_profile(all_dir, profile_fig_dir)
    # analyze_jaccard(fig_dir)