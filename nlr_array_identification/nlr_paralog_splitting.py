#!/usr/bin/env python
# coding: utf-8

# Import packages
import ete3 as ete
from ete3 import PhyloTree
import sys
import numpy as np
import pandas as pd
import os as os

# Print help message
if len(sys.argv)<2:
    print("Usage: python3 paralog_splitting.py tree_file orthogroup_id")
else:
    #if __name__ == "__main__":
    tree_name = sys.argv[1]
    orthogroup = sys.argv[2]

    # Read tree
    t = ete.EvolNode(tree_name)

    # Split tree by duplication events using the SO algorithm
    t.set_species_naming_function(lambda node: node.name.split("_")[0])
    print(t.get_ascii(attributes=["name", "species"], show_internal=False ))

    # Separate new trees and print
    new_trees = []
    i = 0
    original_stdout = sys.stdout
    for node in t.split_by_dups():
        new_trees.append(node)
        i += 1
        og_i = str(orthogroup) + "." + str(i) + "_printTree.txt"
        with open(og_i, mode='wt', encoding='utf-8') as f:
            sys.stdout = f
            print(node)
            sys.stdout = original_stdout

    # Write out new trees as their own orthogroups
    i = 0
    for tree in new_trees:
        i += 1
        og_i = str(orthogroup) + "." + str(i) + ".txt"
        tree_leaves = tree.get_leaf_names()
        with open(og_i, mode='wt', encoding='utf-8') as f:
            f.write('\n'.join(tree_leaves))


