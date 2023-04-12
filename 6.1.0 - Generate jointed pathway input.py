# -*- coding: utf-8 -*-

### Jointed pathway analysis input data generator ###

import os
import pandas as pd

## Folders ##

Folder1 = "Data/Metabolomics/Data analysis"
Folder2 = "Results/Supplementary tables"
Folder3 = "Data/Metabolomics/Jointed pathway analysis data/Metabolyst input"

## Files ##

File1 = "Data analysis overview.xlsx"
File2 = "Supplementary table 1.xlsx"

File3 = "Metabolomics_compound_logFC.txt"
File4 = "RNA_Gene_symbol_logFC.txt"

df1 = pd.read_excel(os.path.join(Folder1,File1))
df2 = pd.read_excel(os.path.join(Folder2,File2), sheet_name="Differentially expressed genes")

# Compound header = #compound	log2FC

df1 = df1[['HMDB','Log2FC']]
df1.rename(columns={'HMDB': "#compound","Log2FC": "log2FC"},inplace=True)
df1 = df1.replace("HMDB0062769","METPA0843")
df1.to_csv(os.path.join(Folder3,File3),index=False,sep="\t")

# Gene header = #Official	log2FC
df2 = df2[['Gene',"log2FoldChange"]]
df2 = df2.round(4)
df2.rename(columns={'Gene': "#Official","log2FoldChange": "log2FC"},inplace=True)
df2.to_csv(os.path.join(Folder3,File4),index=False,sep="\t")
