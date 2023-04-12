# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 09:19:10 2023

@author: dcs839
"""

### Differentially expressed transporters ### 

## Libraries ##

import os
import pandas as pd

## Folders ##

Folder1 = "Results/Supplementary tables"
Folder2 = "Data/RNA/Panther DB/Subcellular location"
Folder3 = "Results/Transport analysis"

os.makedirs(Folder3,exist_ok=True)

## Files ##

File1 = "Supplementary table 1.xlsx"
File2 = "Transporters_with_scores.xlsx"

File3 = "Transport_analysis_overvew_50TPM.xlsx"

cutoff = 50
## Load data ##

df_1 = pd.read_excel(os.path.join(Folder1,File1),sheet_name='Differentially expressed genes')
df_2 = pd.read_excel(os.path.join(Folder2,File2),sheet_name='Transport scores')
df_3 = df_1[df_1['Gene'].isin(df_2['Genes'])].reset_index(drop=True)


df_2_dict = df_2.set_index('Genes')['Plasma membrane scores'].to_dict()
df_3['Plasma membrane scores'] = df_3['Gene'].map(df_2_dict)

#df_3.insert(2, 'Plasma membrane score', df_2['Plasma membrane scores'])


df_up = df_3[(df_3['log2FoldChange'] > 0) & (df_3['Lightphase Mean (TPM)']+df_3['Darkphase Mean (TPM)'] > cutoff)]
df_down = df_3[(df_3['log2FoldChange'] < 0) & (df_3['Lightphase Mean (TPM)']+df_3['Darkphase Mean (TPM)'] > cutoff)]


writer = pd.ExcelWriter(os.path.join(Folder3,File3), engine='xlsxwriter')
df_3.to_excel(writer,sheet_name='All',index=False)
df_up.to_excel(writer,sheet_name='Upregulated',index=False)
df_down.to_excel(writer,sheet_name='Downregulated',index=False)
writer.save()

