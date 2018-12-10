#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 20:27:02 2018

Testing tree formatting

@author: joshuasacher
"""
from Bio import Phylo
from matplotlib import pyplot as plt
import os

t = "tree_all"
tree = Phylo.read(f"{t}_phyml_tree.txt", "newick")
tree = tree.as_phyloxml()
original_id = "NP_008881.2"

# Stylize the tree
# Colorblind-safe colors checked with ColorOracle: http://colororacle.org/
for clade in tree.find_clades():
    # Bold lines
    clade.width = 3
    # Red if known gene or false positive
    if str(clade.name).startswith("gi|"):
        clade.color = "#e4002b"
    # Blue for originally searched gene
    elif clade.name == original_id:
        clade.color = "#006db6"
    # Black for comparitor nodes
    elif clade.name is not None and not clade.color:
        clade.color = "#000000"
    # Gray for non-terminal nodes
    elif not clade.name:
        clade.color = "#63666a"
    # Green for novel genes
    if str(clade.name).endswith("***"):
        clade.color = "#00bf71"

# Plot tree using number of nodes to determine size
tree_len = len(tree.get_terminals())
plt.rc("font", size=18)  # Bigger font for easier reading
fig = plt.figure(figsize=(1.6 * tree_len, tree_len), dpi=300)
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes)
fig.savefig("tree.png", format='png', bbox_inches='tight', dpi=300)
fig.savefig("tree_transparent.png", format='png', bbox_inches='tight', dpi=300,
            transparent=True)
tree_file = t
print(f"Tree images saved as {tree_file + '.png'}"
      f"and {tree_file + '_transparent.png'} to {os.getcwd()}\n")