#!/usr/bin/env python
# coding: utf-8

# run with: python3 py_TEST_LOCAL_commands.py <operational_condition> <MPR>

import sys
import numpy as np
import biotite
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import biotite.sequence.align as align
import biotite.database.entrez as entrez
import pandas as pd
import re as re
import os as os


#if __name__ == "__main__":
filename = sys.argv[1]
outfile = sys.argv[2]
cutoff = int(sys.argv[3])

# Percent ident matrix from clustalw2 -TREE -PIM -BOOTSTRAP=1000 -INFILE=$file
pid_mat = pd.read_table(filename, sep='\t', header = None, index_col = 0)
names = list(pid_mat.index)
pid_mat.columns = names

# Function to get row,col positions of cells over a specific % identity
def getIndexes(dfObj, value):
    ''' Get index positions of value in dataframe i.e. dfObj.'''
    listOfPos = list()
    # Get bool dataframe with True at positions where the given value exists
    result = dfObj > value
    # Get list of columns that contains the value
    seriesObj = result.any()
    columnNames = list(seriesObj[seriesObj == True].index)
    # Iterate over list of columns and fetch the rows indexes where value exists
    for col in columnNames:
        rows = list(result[col][result[col] == True].index)
        for row in rows:
            listOfPos.append((row, col))
    # Return a list of tuples indicating the positions of value in the dataframe
    return listOfPos

# Get list index positions for all values >= cutoff
listOfPositions = getIndexes(pid_mat, cutoff)

#for i in range(len(listOfPositions)):
#    print(listOfPositions[i])

# Turn dictionary of index positions into a dataframe 
paralogs = pd.DataFrame()
paralogs = paralogs.from_dict(listOfPositions)
paralogs.columns = ['GeneA', 'GeneB']

# What genotype does each gene belong to
paralogs['GeneAFam'] = paralogs.GeneA.str.replace('_Chr.*', '', regex=True)
paralogs['GeneBFam'] = paralogs.GeneB.str.replace('_Chr.*', '', regex=True)

# Make sure genes belong to same genotype
paralogs = paralogs[(paralogs['GeneAFam'] == paralogs['GeneBFam']) == True]

# Remove self-self alignments 
paralogs = paralogs[(paralogs['GeneA'] != paralogs['GeneB']) == True]

# Add % identity column to dataframe of index positions
pident = list()
for i in range(len(paralogs)):
    row = paralogs.iloc[i,0]
    col = paralogs.iloc[i,1]
    val = pid_mat.loc[row, col]
    pident.append((val))

paralogs['pident'] = pident
paralogs

# Subset the dataframe, keeping only index positions and % identity
paralogs['pident'] = pident
paralogs = paralogs[["GeneA", "GeneB", "pident"]]
paralogs

# Combine each set of index positions, i.e. row, into list
full_list = tuple(paralogs[['GeneA', 'GeneB']].values.tolist())

# Function to merge lists based on intersections between them
# If one or more element(s) is shared between two  lists, combine them
# Proceed until the lists can no longer be merged
def merge(lsts):
    sets = [set(lst) for lst in lsts if lst]
    merged = True
    while merged:
        merged = False
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = True
                    common |= x
            results.append(common)
        sets = results
    return sets


# Merge lists
merged_full_list = merge(full_list)
merge_full_df = pd.DataFrame(merged_full_list)

# Write output files 
cur_dir = os.getcwd() 
merge_full_df.to_csv(str(cur_dir + "/" + outfile), index = False, header=False)

print(merge_full_df)
