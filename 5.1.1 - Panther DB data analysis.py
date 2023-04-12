# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 13:50:13 2023

@author: dcs839
"""

### Panther DB data analysis ###

## Libraries ##

import os
import pandas as pd

## Folders ##

Folder1 = "Data/RNA/Panther DB/Raw data"
Folder2 = "Data/RNA/Panther DB/Clean data"
Folder5 = "Results/Supplementary tables"
os.makedirs(Folder2,exist_ok=True)

## Files ##

File1 = "Panther_protein_classes.txt"
File2 = "Panther_transporter_list.txt"

File3 = "Enrichment data protein classes.xlsx"
File4 = "Enrichment data Transporters.xlsx"

File5 = "Supplementary table 1.xlsx"

File6 = "Transporter_list_correct.xlsx"

## Variables ##
"""
cutoff = 40

## Load protein class data ##

df_classes = pd.read_csv(os.path.join(Folder1,File1),sep="\t",header=None)
# Filter out unwanted columns #
df_classes = df_classes.drop([0,3,4], axis=1)
# Change names of protein classes (remove parentheses and capitalize) #
df_classes[1] = df_classes[1].map(lambda x: x.split(" (")[0])
df_classes[1] = df_classes[1].map(lambda x: x[0].upper()+x[1:])
# extract not assigned value and change to "Unclassified' #
Not_assigned = df_classes[df_classes[1] == "No PANTHER category is assigned"].reset_index(drop=True)
Not_assigned = Not_assigned.replace({"No PANTHER category is assigned":"Unclassified"})
# collapse all values under cutoff to "Others" #
others = df_classes.loc[df_classes[2] < cutoff,[2]].sum(axis=0).reset_index(drop=True).to_frame()
others.rename(columns={0: 2},inplace=True)
others.insert(0,column=1,value='Others')
# Remove others and unclassified from dataframe, reset index and sort according to decending values #
df_classes = df_classes.drop(df_classes[(df_classes[2] < cutoff) | (df_classes[1] == "No PANTHER category is assigned")].index)
df_classes = df_classes.sort_values(by=[2],ascending=False,ignore_index=True)
# Add others and then unclassified as last values
df_classes = pd.concat([df_classes,others,Not_assigned],ignore_index=True)
# Rename columns #
df_classes.rename(columns={1: "Classes",2:"Values"},inplace=True)

## Save protein classes to excel file ##
writer = pd.ExcelWriter(os.path.join(Folder2,File3), engine='xlsxwriter')
df_classes.to_excel(writer,sheet_name='Enrichment data',index=False)
writer.save()
"""
## load transporter data #
"""
df_transporter = pd.read_csv(os.path.join(Folder1,File2),sep="\t",header=None)
df_transporter = df_transporter.drop([0,2,3,4,5], axis=1)
writer = pd.ExcelWriter(os.path.join(Folder2,File4), engine='xlsxwriter')
df_transporter.to_excel(writer,sheet_name='Transporters',index=False,header=False)
writer.save()
"""
df_transport = pd.read_excel(os.path.join(Folder2,File4),header=None)
df_supp = pd.read_excel(os.path.join(Folder5,File5),sheet_name='Differentially expressed genes')
col_list1 = df_transport[0].values.tolist()
col_list2 = df_supp["Gene"].values.tolist()
df_supp_transport = df_supp[df_supp['Category'] == 'Transporter']
df_supp_channels = df_supp[df_supp['Category'] == 'Channel']
col_list3 = df_supp_transport["Gene"].values.tolist()
col_list4 = df_supp_channels["Gene"].values.tolist()
correct_gene_name = []
for key in col_list1:
    if ',' in key:
        key = key.split(",")
        for item in key:
            if item in col_list2:
                correct_gene_name.append(item)
    else:
        if key in col_list2:
            correct_gene_name.append(key)
        else:
            pass
            #print(key)

addition = list(set(col_list3+col_list4)-set(correct_gene_name))

transporter_list = correct_gene_name+addition

transporter_df = pd.DataFrame({'Genes':transporter_list})
#(list(set(col_list1)-set(col_list3)-set(col_list4)))


transporter_df.to_excel(os.path.join(Folder2,File6),sheet_name="Transporters",index=False)

