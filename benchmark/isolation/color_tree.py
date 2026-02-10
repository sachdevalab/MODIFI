# -*- coding: utf-8 -*-
import subprocess
import os
from pathlib import Path
import sys
import pandas as pd
import re
import seaborn as sns
import matplotlib.pyplot as plt

from sample_object import My_sample, Isolation_sample

BAC_PHYLUM_COLORS = {
    "Pseudomonadota":    "#d8b365",  # tan / brown
    "Bacillota":         "#f46d43",  # orange
    "Bacillota_A":       "#8da0cb",  # blue-lavender
    "Desulfobacterota":       "#fc8d62",  # light orange / salmon (same family as Bacillota)
    "Actinomycetota":    "#66c2a5",  # teal-green
    "Bacteroidota":      "#e78ac3",  # pink
    "Campylobacterota":  "#a6d854",  # light green
    "Acidobacteriota":   "#1b9e77",  # dark green
    "Verrucomicrobiota": "#7570b3",  # muted purple
    "Chloroflexota":   "#e7298a",  # magenta

    "Bacillota_C":    "#b15928",  # brown/rust
    "Nitrospirota":  "#fdbf6f",  # light orange/yellow
    "Bacillota_I":    "#ff7f00",  # bright orange
    "Desulfobacterota_B": "#cab2d6",  # light purple
    "Patescibacteria":     "#6a3d9a",  # deep purple
    "Gemmatimonadota":      "#fb9a99",  # light red/pink
    "Others":            "#e6e6e6"   # very light neutral (only gray)
}



header = """TREE_COLORS
#use this template to define branch colors and styles, colored ranges and label colors/font styles/backgrounds
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
SEPARATOR SPACE
#SEPARATOR COMMA

#First 3 fields define the node, type and color
#Possible types are:
#'range': defines a colored range (colored background for labels/clade)
#'clade': defines color/style for all branches in a clade
#'branch': defines color/style for a single branch
#'label': defines font color/style for the leaf label
#'label_background': defines the leaf label background color

#The following additional fields are required:
#for 'range', field 4 defines the colored range label (used in the legend)

#The following additional fields are optional:
#for 'label', field 4 defines the font style ('normal',''bold', 'italic' or 'bold-italic') and field 5 defines the numeric scale factor for the font size (eg. with value 2, font size for that label will be 2x the standard size)
#for 'clade' and 'branch', field 4 defines the branch style ('normal' or 'dashed') and field 5 defines the branch width scale factor (eg. with value 0.5, branch width for that clade will be 0.5 the standard width)

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR
"""

Label_header = """LABELS
#use this template to change the leaf labels, or define/change the internal node names
#additionally, you can specify a custom class for internal nodes, which can be used to automatically collapse them
#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file.

#SEPARATOR TAB
#SEPARATOR SPACE
SEPARATOR COMMA

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#NODE_ID,LABEL,CLASS

#Examples

#note that the class field is optional

#define a name and class for an internal node. Class 'kingdom' will be available when using the #automatic clade collapsing function
#9031|9606,Metazoa,kingdom

#change the label for leaf node 9606\n"""

gradient_header="""
DATASET_GRADIENT
#In gradient datasets, each ID is associated to a single numeric value which is converted to a colored box based on the gradient defined.

#lines starting with a hash are comments and ignored during parsing

#=================================================================#
#                    MANDATORY SETTINGS                           #
#=================================================================#
#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throught this file.
#SEPARATOR TAB
SEPARATOR SPACE
#SEPARATOR COMMA

#label is used in the legend table (can be changed later)
DATASET_LABEL label 1

#dataset color (can be changed later)
COLOR #ff0000

#=================================================================#
#                    OPTIONAL SETTINGS                            #
#=================================================================#

#=================================================================#
#     all other optional settings can be set or changed later     #
#           in the web interface (under 'Datasets' tab)           #
#=================================================================#


#Each dataset can have a legend, which is defined using LEGEND_XXX fields below
#For each row in the legend, there should be one shape, color and label.
#Optionally, you can define an exact legend position using LEGEND_POSITION_X and LEGEND_POSITION_Y. To use automatic legend positioning, do NOT define these values
#Optionally, shape scaling can be present (LEGEND_SHAPE_SCALES). For each shape, you can define a scaling factor between 0 and 1.
#To order legend entries horizontally instead of vertically, set LEGEND_HORIZONTAL to 1
#Shape should be a number between 1 and 6, or any protein domain shape definition.
#1: square
#2: circle
#3: star
#4: right pointing triangle
#5: left pointing triangle
#6: checkmark


#LEGEND_TITLE Dataset legend
#LEGEND_SCALE 1
#LEGEND_POSITION_X 100
#LEGEND_POSITION_Y 100
#LEGEND_HORIZONTAL 0
#LEGEND_SHAPES 1 2 3
#LEGEND_COLORS #ff0000 #00ff00 #0000ff
#LEGEND_LABELS value1 value2 value3
#LEGEND_SHAPE_SCALES 1 1 0.5

#width of the gradient strip
#STRIP_WIDTH 25

#left margin, used to increase/decrease the spacing to the next dataset. Can be negative, causing datasets to overlap.
#MARGIN 0

#border width; if set above 0, a border of specified width (in pixels) will be drawn around the gradient strip
#BORDER_WIDTH 0

#border color; used when BORDER_WIDTH is above 0
#BORDER_COLOR #0000ff

#automatically create and display a legend based on the color gradient defined below
#AUTO_LEGEND 1

#define the gradient colors. Values in the dataset will be mapped onto the corresponding color gradient.
#COLOR_MIN #ff0000
#COLOR_MAX #0000ff

#you can specify a gradient with three colors (e.g red to yellow to green) by setting 'USE_MID_COLOR' to 1, and specifying the midpoint color
#USE_MID_COLOR 1
#COLOR_MID #ffff00

#always show internal values; if set, values associated to internal nodes will be displayed even if these nodes are not collapsed. It could cause overlapping in the dataset display.
#SHOW_INTERNAL 1

#display or hide the dataset label above the gradient strip
#SHOW_LABELS 1

#text label size factor
#SIZE_FACTOR 1

#text label rotation
#LABEL_ROTATION 0

#text label shift in pixels (positive or negative)
#LABEL_SHIFT 0

#align the dataset label to the tree circle; only applies in circular display mode
#LABEL_ALIGN_TO_TREE,0

#Internal tree nodes can be specified using IDs directly, or using the 'last common ancestor' method described in iTOL help pages
#=================================================================#
#       Actual data follows after the "DATA" keyword              #
#=================================================================#
DATA
#ID2 value2
#9606 10000
#LEAF1|LEAF2 11000"""


def collect_iso_ctgsall_dir(all_dir, iso_genome_list_file, filtered_df):
    iso_genome_list = []
    filtered_df = filtered_df[filtered_df["Motif_Num"] > 0] ## consider only those with motifs

    prefix_list = filtered_df['Sample'].tolist()
    print (f"Collecting isolation genomes for {len(prefix_list)} samples...")
    for prefix in prefix_list:
        isolation_obj = Isolation_sample(prefix, all_dir)
        isolation_obj.read_depth()
        isolation_obj.read_MGE()
        genome_list, contig_list = isolation_obj.get_iso_good_ctgs(min_depth=10, min_len=500000)
        iso_genome_list += genome_list
    with open(iso_genome_list_file, "w") as f:
        for genome in iso_genome_list:
            f.write(genome + "\n")
    print (f"Total isolation contigs collected: {len(iso_genome_list)}")
  
def single_run(resultdir, genome_dir):
    data = []
    archea_list = ["SRR27457941", "SRR31014709"]

    for folder in os.listdir(resultdir):
        prefix = folder
        if prefix in archea_list:
            continue
        sample_obj = Isolation_sample(prefix, resultdir)

        if not os.path.exists(sample_obj.gtdb):
            print (f"GTDB file not found for {prefix}, skipping...")
            continue
        pure_anno = sample_obj.check_pure2()

        sample_obj.get_phylum()
        motif_num, unique_motifs = sample_obj.get_unique_motifs(min_frac=0.3, min_sites = 100)
        mge_bool = sample_obj.get_MGE_bool()
        average_dp = sample_obj.get_average_depth()

        data.append([prefix, sample_obj.lineage, motif_num, mge_bool, pure_anno, average_dp, sample_obj.reference_fasta])

    df = pd.DataFrame(data, columns=["Sample", "Lineage", "Motif_Num", "MGE_bool", "Pure_anno", "Average_DP", "genome"])
    df.to_csv(f"{tree_results}/isolation_sample_summary.tsv", sep="\t", index=False)
    df = filter_df(df, min_dp = 10)
    run_taxa_dict = {}
    sample_meta_dict = {}
    for index, row in df.iterrows():
        run_taxa_dict[row["Sample"]] = row["Lineage"]
        sample_meta_dict[row["Sample"]] = {
            "motif_num": row["Motif_Num"],
            "mge_bool": row["MGE_bool"],
            "average_dp": row["Average_DP"]
        }
        ## ln -s genome to genome_dir
        src = row["genome"]
        dst = os.path.join(genome_dir, os.path.basename(row["Sample"] + ".fa"))
        if not os.path.exists(dst):
            os.symlink(src, dst)
        
    return run_taxa_dict, sample_meta_dict, df

def filter_df(df, min_dp = 10):
    """
    count how many in original df, and filter by pure, and filter by min_dp
    And after pure, how many have motif_num > 0, and how many have MGE_bool = True 
    """
    original_num = df.shape[0]
    print (f"Original samples: {original_num}")
    df_pure = df[df['Pure_anno'] == "pure"]
    pure_num = df_pure.shape[0]
    print (f"Pure samples: {pure_num} ({pure_num/original_num:.2%})")
    df_dp = df_pure[df_pure['Average_DP'] >= min_dp]
    dp_num = df_dp.shape[0]
    print (f"Pure samples with DP >= {min_dp}: {dp_num} ({dp_num/pure_num:.2%})")
    motif_positive_num = df_dp[df_dp['Motif_Num'] > 0].shape[0]
    print (f"Pure samples with DP >= {min_dp} and Motif_Num > 0: {motif_positive_num} ({motif_positive_num/dp_num:.2%})")
    mge_positive_num = df_dp[df_dp['MGE_bool'] == 1].shape[0]
    print (f"Pure samples with DP >= {min_dp} and MGE_bool = 1: {mge_positive_num} ({mge_positive_num/dp_num:.2%})")
    ## count how many unique species and phylum in df_dp
    unique_species = set()
    unique_phylum = set()
    no_species_count = 0
    for index, row in df_dp.iterrows():
        lineage = row["Lineage"]
        species = lineage.split(";")[-1][3:] if ";" in lineage else "Unclassified"
        phylum = lineage.split(";")[1][3:] if ";" in lineage else "Unclassified"
        if species == "":
            # print (lineage)
            # print (row)
            species = "Unclassified"
            no_species_count += 1
        else:
            unique_species.add(species)
        unique_phylum.add(phylum)
    # print (unique_species)
    print (f"Unique species in filtered samples: {len(unique_species)}")
    print (f"Unique phylum in filtered samples: {len(unique_phylum)}")
    print (f"Samples without species information: {no_species_count}")

    ## count the proportion of samples in each phylum with motif_num > 0
    print ("\nProportion of samples with Motif_Num > 0 in each phylum:")
    phylum_groups = df_dp.groupby(df_dp['Lineage'].apply(lambda x: x.split(";")[1][3:] if ";" in x else "Unclassified"))
    for phylum, group in phylum_groups:
        total_count = group.shape[0]
        motif_positive_count = group[group['Motif_Num'] > 0].shape[0]
        proportion = motif_positive_count / total_count if total_count > 0 else 0
        print (f"{phylum}: {motif_positive_count}/{total_count} ({proportion:.2%})")
    
    return df_dp

def color_phylum(run_taxa_dict, tree_results, sample_meta_dict):
    ## in iTol, color by phylum using general colors
    ## collect phyllums in run_taxa_dict, and assign colors
    taxa_list = set(list(run_taxa_dict.values()))
    
    # Extract unique phylums
    phylums = set()
    for taxa in run_taxa_dict.values():
        if taxa != "unknown":
            phylum = taxa.split(";")[1] if len(taxa.split(";")) > 1 else "unknown"
            phylum = phylum.replace("p__", "")
            phylums.add(phylum)
    
    # Use predefined PHYLUM_COLORS
    phylum_list = sorted(list(phylums))  # Sort for consistent coloring
    
    # Assign colors from PHYLUM_COLORS dictionary
    phylum_color = {}
    for phylum in phylum_list:
        if phylum in BAC_PHYLUM_COLORS:
            phylum_color[phylum] = BAC_PHYLUM_COLORS[phylum]
        else:
            # Use "Others" color for phylums not in the dictionary
            phylum_color[phylum] = BAC_PHYLUM_COLORS["Others"]
    
    # Add default color for unknown (use same as Others)
    phylum_color["unknown"] = BAC_PHYLUM_COLORS["Others"]
    
    print(f"Found {len(phylum_list)} unique phylums using predefined colors:")
    for phylum, color in phylum_color.items():
        print(f"  {phylum}: {color}")
    Label_anno = Label_header
    

    motif_num_anno = gradient_header + "\n"
    mge_num_anno = gradient_header + "\n"
    depth_anno = gradient_header + "\n"
    with open(f"{tree_results}/phylum_annotations.txt", "w") as f:
        f.write(f"{header}\n")
        for bin_name, taxa in run_taxa_dict.items():
            phylum = taxa.split(";")[1] if len(taxa.split(";")) > 1 else "unknown"
            phylum = phylum.replace("p__", "")
            color = phylum_color.get(phylum, "#cccccc")
            
            # Extract species information (typically at index 6 in GTDB taxonomy)
            species = taxa.split(";")[6] if len(taxa.split(";")) > 6 else "unknown"
            species = species.replace("s__", "")

            f.write(f"{bin_name} range {color} {phylum}\n")
            node_rename = f"{species}_{bin_name}"
            Label_anno += f"{bin_name},{node_rename}\n"
            if sample_meta_dict[bin_name]['motif_num'] is not None:
                try:
                    motif_num = float(sample_meta_dict[bin_name]['motif_num'])
                    if motif_num > 1:
                        motif_num = 2
                    motif_num_anno += f"{bin_name} {motif_num}\n"
                except (ValueError, TypeError):
                    print(f"Warning: could not convert motif_num to float for {bin_name}: {sample_meta_dict[bin_name]['motif_num']}")
                    motif_num_anno += f"{bin_name} 0\n"
            
            if sample_meta_dict[bin_name]['mge_bool'] is not None:
                try:
                    mge_bool = float(sample_meta_dict[bin_name]['mge_bool'])
                    mge_num_anno += f"{bin_name} {mge_bool}\n"
                except (ValueError, TypeError):
                    print(f"Warning: could not convert mge_bool to float for {bin_name}: {sample_meta_dict[bin_name]['mge_bool']}")
                    mge_num_anno += f"{bin_name} 0\n"
            
            try:
                average_dp = float(sample_meta_dict[bin_name]['average_dp'])
                if average_dp >= 10:
                    depth_anno += f"{bin_name} 1\n"
                else:
                    depth_anno += f"{bin_name} 0\n"
            except (ValueError, TypeError):
                print(f"Warning: could not convert average_dp to float for {bin_name}: {sample_meta_dict[bin_name]['average_dp']}")
                depth_anno += f"{bin_name} 0\n"
            # Add to legend data
            # legend_data += f"{bin_name} {color} {phylum}\n"
    

    print(f"\nGenerated iTOL annotation file with {len(run_taxa_dict)} entries using Set1 colors")
    print(f"Legend includes {len([p for p in phylum_color.keys() if p != 'unknown'])} phylums")
    # print(run_taxa_dict)
    
    ## save label annotation file
    with open(f"{tree_results}/label_annotations.csv", "w") as f:
        f.write(Label_anno)
    with open(f"{tree_results}/motif_num_annotations.txt", "w") as f:
        f.write(motif_num_anno)
    with open(f"{tree_results}/mge_num_annotations.txt", "w") as f:
        f.write(mge_num_anno)
    with open(f"{tree_results}/depth_annotations.txt", "w") as f:
        f.write(depth_anno)


if __name__ == "__main__":
    resultdir = f"/home/shuaiw/borg/paper/isolation//batch2_results/"
    tree_results = "/home/shuaiw/borg/paper/isolation//GTDB_tree/anno/"
    genome_dir = "/home/shuaiw/borg/paper/isolation//GTDB_tree/genomes"
    genome_list = "/home/shuaiw/borg/paper/specificity/iso_genome.list"  ## for drep
    run_taxa_dict, sample_meta_dict, filtered_df = single_run(resultdir, genome_dir) ## collect isolation genomes with high dp
    # collect_iso_ctgsall_dir(resultdir, genome_list, filtered_df)
    # color_phylum(run_taxa_dict, tree_results, sample_meta_dict)