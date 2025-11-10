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

  
def single_run(resultdir):
    
    run_taxa_dict = {}
    sample_meta_dict = {}
    for folder in os.listdir(resultdir):
        prefix = folder
        sample_obj = Isolation_sample(prefix, resultdir)

        if not os.path.exists(sample_obj.gtdb):
            print (f"GTDB file not found for {prefix}, skipping...")
            continue
        sample_obj.get_phylum()
        run_taxa_dict[prefix] = sample_obj.lineage
        motif_num, unique_motifs = sample_obj.get_unique_motifs()
        sample_meta_dict[prefix] = {
            "motif_num": motif_num
        }
    print ("No. of runs (samples) processed:", len(run_taxa_dict))
    return run_taxa_dict, sample_meta_dict

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
    
    # Use seaborn Set1 color palette
    phylum_list = sorted(list(phylums))  # Sort for consistent coloring
    n_colors = len(phylum_list)
    
    # Get Set1 palette colors (Set1 has 9 distinct colors)
    set1_colors = list(sns.color_palette("Set2", n_colors=max(n_colors, 9)))
    
    # Convert RGB to hex format for iTOL
    phylum_color = {}
    for i, phylum in enumerate(phylum_list):
        if i < len(set1_colors):
            rgb = set1_colors[i]
            hex_color = "#{:02x}{:02x}{:02x}".format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            phylum_color[phylum] = hex_color
        else:
            # If more phylums than Set1 colors, cycle through
            rgb = set1_colors[i % len(set1_colors)]
            hex_color = "#{:02x}{:02x}{:02x}".format(int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            phylum_color[phylum] = hex_color
    
    # Add default color for unknown
    phylum_color["unknown"] = "#cccccc"
    
    print(f"Found {len(phylum_list)} unique phylums using Seaborn Set1 palette:")
    for phylum, color in phylum_color.items():
        print(f"  {phylum}: {color}")
    Label_anno = Label_header
    

    motif_num_anno = gradient_header + "\n"
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
                motif_num_anno += f"{bin_name} {sample_meta_dict[bin_name]['motif_num']}\n"
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


if __name__ == "__main__":
    resultdir = f"/groups/banfield/projects/multienv/methylation_temp/batch2_results/"
    tree_results = "/groups/banfield/projects/multienv/methylation_temp/GTDB_tree/anno/"
    run_taxa_dict, sample_meta_dict = single_run(resultdir)
    color_phylum(run_taxa_dict, tree_results, sample_meta_dict)