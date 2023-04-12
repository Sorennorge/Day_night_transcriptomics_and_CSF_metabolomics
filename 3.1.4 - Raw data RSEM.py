# -*- coding: utf-8 -*-

### Convert all RSEM output files to csv, extracting TPM column ###

## Libraries ##

import os

## Folders ##

Folder_in = "Data/RNA/Raw data/RSEM"
Folder_out = "Data/RNA/Raw data/TPM"

os.makedirs(Folder_out, exist_ok=True)

### Prep raw data to standard input ###

for x in range(1,13,1):
    Ensembl_counts = {}
    with open(os.path.join(Folder_in,"RSEM_sample_{}.txt".format(x)),'r') as read:
        next(read)
        for line in read:
            line = line.strip().split("\t")
            Ensembl_counts[line[0]] = line[5]
    read.close
    
    ## Save data as ensembl plus TPM
    with open(os.path.join(Folder_out,"Sample_{}_TPM.csv".format(x)),'w+') as out:
        out.write("Entry;TPM\n")
        for key in Ensembl_counts:
            out.write("{};{}\n".format(key,Ensembl_counts[key]))
    out.close
