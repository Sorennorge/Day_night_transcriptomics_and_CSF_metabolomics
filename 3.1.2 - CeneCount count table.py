# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 08:22:17 2021

@author: dcs839
"""

### Generate count table from gene count files (RNA-STAR) ###

## Import libraries ##

import os

## Folders ##

Folder = "Data/RNA//Raw data/Raw_counts"
Folder_out = "Data/RNA/Count Tables"

os.makedirs(Folder_out, exist_ok=True)

## Files ##

File_out = "Count_table_Sleep_RNASTAR.csv"

## Variables

Gene_table = {}
    
Sleep_header = []
for x in range(1,7,1):
    Sleep_header.append("Day {}".format(x))
for x in range(1,7,1):
    Sleep_header.append("Night {}".format(x))

## create count table ###

for x in range(1,13,1):
    # Go through samples 1-12
    Gene_counts = {}
    with open(os.path.join(Folder,"Sample_{}_rawcounts.csv".format(x)),'r') as read:
        next(read)
        for line in read:
            line = line.strip().split(";")
            if line[0] not in Gene_table:
                #if ensembl not in gene table create array of 
                Gene_table[line[0]] = [0] * 12
                Gene_table[line[0]][x-1] += int(line[1])
            else:
                Gene_table[line[0]][x-1] += int(line[1])
    read.close

## Save count table to file ##
with open(os.path.join(Folder_out,File_out),'w+') as out:
    out.write("Gene;{}\n".format(";".join(Sleep_header)))
    for key in sorted(Gene_table):
        out.write("{};{}\n".format(key,";".join(map(str,Gene_table[key]))))
out.close
