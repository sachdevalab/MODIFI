#!/usr/bin/env python3
"""
Generate a fake phylogenetic tree for iTOL with phyla as leaves.
Each phylum from the PHYLUM_COLORS dictionary will be a leaf in the tree.
"""

import random


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

def generate_balanced_tree(taxa_list, branch_length=1.0):
    """
    Generate a balanced binary tree in Newick format.
    
    Args:
        taxa_list: List of taxon names (leaves)
        branch_length: Branch length for the tree
    
    Returns:
        Newick format tree string
    """
    if len(taxa_list) == 1:
        return f"{taxa_list[0]}:{branch_length}"
    
    # Split taxa into two groups
    mid = len(taxa_list) // 2
    left = taxa_list[:mid]
    right = taxa_list[mid:]
    
    # Recursively build subtrees
    left_tree = generate_balanced_tree(left, branch_length)
    right_tree = generate_balanced_tree(right, branch_length)
    
    return f"({left_tree},{right_tree}):{branch_length}"


def generate_random_tree(taxa_list, min_branch=0.1, max_branch=1.0):
    """
    Generate a random tree in Newick format with variable branch lengths.
    
    Args:
        taxa_list: List of taxon names (leaves)
        min_branch: Minimum branch length
        max_branch: Maximum branch length
    
    Returns:
        Newick format tree string
    """
    if len(taxa_list) == 1:
        bl = random.uniform(min_branch, max_branch)
        return f"{taxa_list[0]}:{bl:.4f}"
    
    if len(taxa_list) == 2:
        bl1 = random.uniform(min_branch, max_branch)
        bl2 = random.uniform(min_branch, max_branch)
        return f"({taxa_list[0]}:{bl1:.4f},{taxa_list[1]}:{bl2:.4f})"
    
    # Randomly split taxa
    random.shuffle(taxa_list)
    split_point = random.randint(1, len(taxa_list) - 1)
    left = taxa_list[:split_point]
    right = taxa_list[split_point:]
    
    # Recursively build subtrees
    left_tree = generate_random_tree(left, min_branch, max_branch)
    right_tree = generate_random_tree(right, min_branch, max_branch)
    
    bl = random.uniform(min_branch, max_branch)
    return f"({left_tree},{right_tree}):{bl:.4f}"


def generate_color_annotation(phyla_colors):
    """
    Generate iTOL color annotation file for the phyla.
    
    Args:
        phyla_colors: Dictionary mapping phylum names to hex colors
    
    Returns:
        Color annotation string for iTOL
    """
    header = """TREE_COLORS
SEPARATOR TAB
DATA
"""
    lines = [header]
    
    # Add color for each leaf (phylum)
    for phylum, color in phyla_colors.items():
        # Label color
        lines.append(f"{phylum}\tlabel\t{color}\tnormal\t1")
        # Branch color (clade with single leaf)
        lines.append(f"{phylum}\tbranch\t{color}\tnormal\t2")
    
    return "\n".join(lines)


def generate_label_annotation(phyla_list):
    """
    Generate iTOL label annotation file.
    
    Args:
        phyla_list: List of phylum names
    
    Returns:
        Label annotation string for iTOL
    """
    header = """LABELS
SEPARATOR TAB
DATA
"""
    lines = [header]
    
    for phylum in phyla_list:
        lines.append(f"{phylum}\t{phylum}")
    
    return "\n".join(lines)


def generate_colored_ranges(phyla_colors):
    """
    Generate iTOL colored ranges annotation.
    
    Args:
        phyla_colors: Dictionary mapping phylum names to hex colors
    
    Returns:
        Colored ranges annotation string for iTOL
    """
    header = """DATASET_COLORSTRIP
SEPARATOR TAB
DATASET_LABEL\tPhylum
COLOR\t#ff0000

LEGEND_TITLE\tPhylum
LEGEND_SHAPES\t"""
    
    shapes = "\t".join(["1"] * len(phyla_colors))
    header += shapes + "\n"
    header += "LEGEND_COLORS\t" + "\t".join(phyla_colors.values()) + "\n"
    header += "LEGEND_LABELS\t" + "\t".join(phyla_colors.keys()) + "\n"
    header += "\nDATA\n"
    
    lines = [header]
    
    for phylum, color in phyla_colors.items():
        lines.append(f"{phylum}\t{color}\t{phylum}")
    
    return "\n".join(lines)


def main():
    import os
    
    # Set output directory
    output_dir = "/home/shuaiw/borg/paper/specificity/fake_tree/"
    os.makedirs(output_dir, exist_ok=True)
    
    # Set random seed for reproducibility
    random.seed(42)
    PHYLA = list(BAC_PHYLUM_COLORS .keys())
    
    # Generate balanced tree
    phyla_sorted = sorted(PHYLA)
    balanced_tree = generate_balanced_tree(phyla_sorted, branch_length=1.0)
    balanced_newick = balanced_tree + ";"
    
    # Generate random tree
    phyla_copy = PHYLA.copy()
    random_tree = generate_random_tree(phyla_copy, min_branch=0.5, max_branch=2.0)
    random_newick = random_tree + ";"
    
    # Write balanced tree
    balanced_path = os.path.join(output_dir, "phylum_tree_balanced.nwk")
    with open(balanced_path, "w") as f:
        f.write(balanced_newick + "\n")
    print(f"Generated: {balanced_path}")
    
    # Write random tree
    random_path = os.path.join(output_dir, "phylum_tree_random.nwk")
    with open(random_path, "w") as f:
        f.write(random_newick + "\n")
    print(f"Generated: {random_path}")
    
    # Generate color annotation
    color_annotation = generate_color_annotation(BAC_PHYLUM_COLORS )
    color_path = os.path.join(output_dir, "phylum_colors.txt")
    with open(color_path, "w") as f:
        f.write(color_annotation)
    print(f"Generated: {color_path}")
    
    # Generate label annotation
    label_annotation = generate_label_annotation(PHYLA)
    label_path = os.path.join(output_dir, "phylum_labels.txt")
    with open(label_path, "w") as f:
        f.write(label_annotation)
    print(f"Generated: {label_path}")
    
    # Generate colored ranges
    colored_ranges = generate_colored_ranges(BAC_PHYLUM_COLORS )
    colorstrip_path = os.path.join(output_dir, "phylum_colorstrip.txt")
    with open(colorstrip_path, "w") as f:
        f.write(colored_ranges)
    print(f"Generated: {colorstrip_path}")
    
    print("\nFiles generated successfully!")
    print("\nTo use in iTOL:")
    print("1. Upload phylum_tree_balanced.nwk or phylum_tree_random.nwk to iTOL")
    print("2. Drag and drop phylum_colors.txt to add branch/label colors")
    print("3. Drag and drop phylum_colorstrip.txt to add colored strip annotation")


if __name__ == "__main__":
    main()
