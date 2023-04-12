# -*- coding: utf-8 -*-

### Evaluate DE analysis of Lightphase/Darkphase rats -- Add information and create supplementary tables ###

### Import libraries ###

import os
import pandas as pd
import numpy as np

## Folders ##

Folder1 = "Data/RNA/Gene info/Biomart"
Folder2 = "Results/RNA/Differential expressed tables"
Folder3 = "Data/RNA/Gene info/Categories"
Folder4 = "Data/RNA/Count Tables"

Folder5 = "Results/Supplementary tables"

os.makedirs(Folder5,exist_ok=True)

## Files ##

# file gathered from https://may2021.archive.ensembl.org/biomart/martview/bf2e14cdc13a5a35276820e26937a5a2
# with only ensembl id and gene symbol / name
Biomart = "BioMart_Rnor6.0.txt"

# Differential expression output files #
Deseq_DE_file = "DEseq2_DE_genes_Analysis_Sleep.csv"
Deseq_all_file = "DEseq2_all_genes_Analysis_Sleep.csv"

# Category file #
# Category file are from earlier work: https://fluidsbarrierscns.biomedcentral.com/articles/10.1186/s12987-022-00335-x
Category_file = "Supplementary Tables.xlsx"

# Count table #

Count_table_file = "Count_table_Sleep_RSEM.csv"

# Output file #
Output_file = "Supplementary table 1.xlsx"

## Variables ##

Ensembl_to_genes = {}

Categories_dict = {}

TPM_dict_Lightphase_mean = {}
TPM_dict_Lightphase_std = {}
TPM_dict_Darkphase_mean = {}
TPM_dict_Darkphase_std = {}

Ensembl_to_category = {}

### Read files ###

## Biomart ##
print("Reading infomation..")
with open(os.path.join(Folder1,Biomart),'r') as read:
    next(read)
    for line in read:
        line = line.strip().split(",")
        if line[0] not in Ensembl_to_genes:
            if line[1] == '':
                Ensembl_to_genes[line[0]] = 'N/A'
            else:
                Ensembl_to_genes[line[0]] = line[1]
read.close

## DE file ##

deseq_df_de = pd.read_csv(os.path.join(Folder2,Deseq_DE_file),
                          sep=";",
                          decimal=",",
                          usecols = [0,1,2,5,6]).rename(columns={'Unnamed: 0': "Ensembl ID"}).set_index('Ensembl ID')
DE_Ensembl_ids = deseq_df_de.index.values.tolist()

## ALL file ##

deseq_df_all = pd.read_csv(os.path.join(Folder2,Deseq_all_file),
                          sep=";",
                          decimal=",",
                          usecols = [0,1,2,5,6]).rename(columns={'Unnamed: 0': "Ensembl ID"}).set_index('Ensembl ID')
ALL_Ensembl_ids = deseq_df_all.index.values.tolist()
print("Done.")

## Category files ##
print("Reading categories..")
# Transport #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S1 - Transporters-pumps',skiprows=[0,1])
Categories_dict['Transporter'] = list(df['Ensembl ID'])
# Channels #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S2 - Channels',skiprows=[0,1])
Categories_dict['Channel'] = list(df['Ensembl ID'])
# GPCRs #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S8 - GPCRs',skiprows=[0,1])
Categories_dict['GPCR'] = list(df['Ensembl ID'])
# RTKs #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S9 - RTKs',skiprows=[0,1])
Categories_dict['RTK'] = list(df['Ensembl ID'])
# Kinases #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S10 - Kinases',skiprows=[0,1])
Categories_dict['Kinase'] = list(df['Ensembl ID'])
# Phosphatases #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S11 - Phosphatases',skiprows=[0,1])
Categories_dict['Phosphatase'] = list(df['Ensembl ID'])
# PDEs #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S12 - PDEs',skiprows=[0,1])
Categories_dict['PDE'] = list(df['Ensembl ID'])
# Cyclases #
df = pd.read_excel(os.path.join(Folder3,Category_file), sheet_name='Table S13 - Cyclases',skiprows=[0,1])
Categories_dict['Cyclase'] = list(df['Ensembl ID'])
print("Assigning categories..")
# Collect ensembl ids and assign categories #
Ensembl_to_category = {}
for key_id in ALL_Ensembl_ids:
    for key in Categories_dict:
        if key_id in Categories_dict[key]:
            if key_id in Ensembl_to_category:
                print(key_id)
                break
            else:
                Ensembl_to_category[key_id] = key
        else:
            pass
print("Done.")      
print("Reading count tables..")  
# Find TPM mean & SD
with open(os.path.join(Folder4,Count_table_file),'r') as read:
    next(read)
    for line in read:
        line = line.strip().split(";")
        # Init #
        Key = line[0]
        Lightphase_array = np.array(line[1:7],dtype=float)
        Darkphase_array = np.array(line[7:13],dtype=float)
        # Lightphase #
        Lightphase_mean = np.round(np.mean(Lightphase_array),2)
        Lightphase_std = np.round(np.std(Lightphase_array, ddof=1),2)
        TPM_dict_Lightphase_mean[line[0]] = Lightphase_mean
        TPM_dict_Lightphase_std[line[0]] = Lightphase_std
        # Darkphase #
        Darkphase_mean = np.round(np.mean(Darkphase_array),2)
        Darkphase_std = np.round(np.std(Darkphase_array, ddof=1),2)
        TPM_dict_Darkphase_mean[line[0]] = Darkphase_mean
        TPM_dict_Darkphase_std[line[0]] = Darkphase_std
        
read.close
print("Done.")

print("Creating dataframe..")
### Create dataframes for output files ###
# Create data frames #
df_de_genes = pd.DataFrame(DE_Ensembl_ids, columns=['Ensembl ID'])
df_all_genes = pd.DataFrame(ALL_Ensembl_ids, columns=['Ensembl ID'])
# Add gene names #
df_de_genes['Gene'] = df_de_genes['Ensembl ID'].map(Ensembl_to_genes)
df_all_genes['Gene'] = df_all_genes['Ensembl ID'].map(Ensembl_to_genes)
## Add Meanbase, Log2FC, Pvalue, and Padj ##
df_de_genes = df_de_genes.join(deseq_df_de, on='Ensembl ID')
df_all_genes = df_all_genes.join(deseq_df_all, on='Ensembl ID')
## Add Category ##
df_de_genes['Category'] = df_de_genes['Ensembl ID'].map(Ensembl_to_category)
df_all_genes['Category'] = df_all_genes['Ensembl ID'].map(Ensembl_to_category)
## Add TPM mean & std for Lightphase and Darkphase ##
# Differentially expressed genes #
df_de_genes['Lightphase Mean (TPM)'] = df_de_genes['Ensembl ID'].map(TPM_dict_Lightphase_mean)
df_de_genes['Lightphase SD (TPM)'] = df_de_genes['Ensembl ID'].map(TPM_dict_Lightphase_std)
df_de_genes['Darkphase Mean (TPM)'] = df_de_genes['Ensembl ID'].map(TPM_dict_Darkphase_mean)
df_de_genes['Darkphase SD (TPM)'] = df_de_genes['Ensembl ID'].map(TPM_dict_Darkphase_std)
# all genes #
df_all_genes['Lightphase Mean (TPM)'] = df_all_genes['Ensembl ID'].map(TPM_dict_Lightphase_mean)
df_all_genes['Lightphase SD (TPM)'] = df_all_genes['Ensembl ID'].map(TPM_dict_Lightphase_std)
df_all_genes['Darkphase Mean (TPM)'] = df_all_genes['Ensembl ID'].map(TPM_dict_Darkphase_mean)
df_all_genes['Darkphase SD (TPM)'] = df_all_genes['Ensembl ID'].map(TPM_dict_Darkphase_std)
print("Done.")

print("Saving dataframes as supplementary table")
### Save dataframes to supplementary file ###
writer = pd.ExcelWriter(os.path.join(Folder5,Output_file), engine='xlsxwriter')
df_all_genes.to_excel(writer,sheet_name='All genes',index=False)
df_de_genes.to_excel(writer,sheet_name='Differentially expressed genes',index=False)
writer.save()
print("All done.")
