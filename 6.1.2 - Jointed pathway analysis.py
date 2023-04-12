# -*- coding: utf-8 -*-

### Jointed pathway analysis - Fold enrichement ###

## Libraries ##
import os
import pandas as pd

## Folders ##

Folder1 = "Data/Metabolomics/Jointed pathway analysis data/Clean data"
Folder2 = "Data/Metabolomics/Pathway data"
Folder3 = "Results/Jointed pathway"
os.makedirs(Folder3,exist_ok=True)

## Files ##

File1 = "Jointed_Pathway_data.xlsx"
File2 = "Pathway reulsts.xlsx"
File3 = "Pathway data.xlsx"

File4 = "Pathway_analysis_overview(Metabolites_plus_jointed).xlsx"

## Load data ##

df_data = pd.read_excel(os.path.join(Folder1,File1),dtype=str)
df_data = df_data[df_data['Reoccurrence'] == 'Yes']

df_results = pd.read_excel(os.path.join(Folder1,File2))
df_results.rename(columns={"Pathways": "Pathway"},inplace=True)
df_results = df_results[df_results['Pathway'].isin(df_data['Pathways'])]

df_pathway = pd.read_excel(os.path.join(Folder2,File3))
df_pathway.rename(columns={"Count": "Hits","Expected (metaboanalyst)": "Expected","Raw P":"Raw p"},inplace=True)
df_pathway = df_pathway.drop(["Ratio metab",], axis=1)

df_jointed_overview = df_results[['Pathway','Impact','Raw p','FDR','Total','Hits','Expected']]
df_pathway_overview = df_pathway[['Pathway','Impact','Raw p','FDR','Total','Hits','Expected','Fold enrichment']]
df_jointed_overview['Fold enrichment'] = df_jointed_overview.apply(lambda row: round(row.Hits/row.Expected,4), axis=1)

df_all = pd.concat([df_pathway_overview.set_index('Pathway'),df_jointed_overview.set_index('Pathway')], axis=1).reset_index()

writer = pd.ExcelWriter(os.path.join(Folder3,File4), engine='xlsxwriter')
df_all.to_excel(writer,sheet_name='Pathway analysis overview',index=False)
writer.save()
