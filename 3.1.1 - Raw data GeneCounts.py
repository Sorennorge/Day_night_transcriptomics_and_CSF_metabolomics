# -*- coding: utf-8 -*-

### Convert RNA-STAR gene count files to csv ###

## Import libraries ##

import os

## Folders #
Folder_in = "Data/RNA/Raw data/GeneCounts"
Folder_out = "Data/RNA/Raw data/Raw_counts"

os.makedirs(Folder_out, exist_ok=True)

### Prep raw data to standard input ###

for x in range(1,13,1):
    Ensembl_counts = {}
    with open(os.path.join(Folder_in,"GeneCount_Sample_{}.txt".format(x)),'r') as read:
        for _ in range(0,4,1):
            next(read)
        for line in read:
            line = line.strip().split("\t")
            Ensembl_counts[line[0]] = int(line[2])
    read.close
    ## Save data as ensembl plus raw data counts
    with open(os.path.join(Folder_out,"Sample_{}_rawcounts.csv".format(x)),'w+') as out:
        out.write("Entry;raw_count\n")
        for key in Ensembl_counts:
            out.write("{};{}\n".format(key,Ensembl_counts[key]))
    out.close

