# -*- coding: utf-8 -*-

### Jointed pathway analysis ###

import os
import pandas as pd

## Folders ##

Folder1 = "Data/Metabolomics/Jointed pathway analysis data/Raw data"
Folder2 = "Data/Metabolomics/Jointed pathway analysis data/Clean data"
os.makedirs(Folder2,exist_ok=True)
Folder3 = "Data/Metabolomics/Pathway data"

## Files ##

File1 = "name_map.csv"
File2 = "gene_name_map.csv"
File3 = "jointpa_matched_features.csv"
File4 = "Pathway data.xlsx"
File5 = "MetaboAnalyst_result_pathway.csv"

File1_clean = "Metabolite_name_mapping.xlsx"
File2_clean = "Gene_name_mapping.xlsx"
File3_clean = "Jointed_Pathway_data.xlsx"
File4_clean = "Pathway reulsts.xlsx"


## Load raw data -> save as clean data ##

# Load metabolite name mapping #
df_metabolite = pd.read_csv(os.path.join(Folder1,File1))
df_metabolite = df_metabolite.drop(["HMDB","PubChem","SMILES","Comment"], axis=1)
df_metabolite = df_metabolite.dropna()
df_metabolite.rename(columns={"Query": "HMDB","Match":"Name","KEGG": "Compound"},inplace=True)
df_metabolite = df_metabolite[["Compound", "HMDB", "Name"]]
# Save clean metabolite data #
writer1 = pd.ExcelWriter(os.path.join(Folder2,File1_clean), engine='xlsxwriter')
df_metabolite.to_excel(writer1,sheet_name='Metabolites',index=False)
writer1.save()

# Load gene name mapping #
df_gene = pd.read_csv(os.path.join(Folder1,File2),dtype=str)
df_gene = df_gene.drop(["Symbol","Name","Comment"], axis=1)
df_gene = df_gene.dropna()
df_gene.rename(columns={"Query": "Gene"},inplace=True)
df_gene = df_gene[["Entrez", "Gene"]]
# Save clean gene mapping data #
writer2 = pd.ExcelWriter(os.path.join(Folder2,File2_clean), engine='xlsxwriter')
df_gene.to_excel(writer2,sheet_name='Genes',index=False)
writer2.save()

df_pathways = pd.read_csv(os.path.join(Folder1,File3))
df_pathways.rename(columns={"Unnamed: 0": "Pathways","matched_features":"Features"},inplace=True)
df_pathways['Features'] = df_pathways.Features.apply(lambda x: x.replace(" ","").split(';'))
df_pathways = df_pathways.replace({"Glycolysis or Gluconeogenesis":"Glycolysis / Gluconeogenesis"})
pathway_dict = df_pathways.set_index('Pathways')['Features'].to_dict()
Gene_features = {}
Compound_features = {}

for key in pathway_dict:
    Gene_features[key] = []
    Compound_features[key] = []
    for item in pathway_dict[key]:
        item = item.split(":")
        if item[0] == 'rno':
            Gene_features[key].append(item[1])
        elif item[0] == 'cpd':
            Compound_features[key].append(item[1])
        else:
            print("Error: incorrect gene or compund tag in list for key: {}, {}".format(key,item))
for key in pathway_dict:
    pathway_dict[key] = "|".join(pathway_dict[key])

for key in Gene_features:
    Gene_features[key] = "|".join(Gene_features[key])
for key in Compound_features:
    Compound_features[key] = "|".join(Compound_features[key])

df_pathway_metabolite_only = pd.read_excel(os.path.join(Folder3,File4), sheet_name="Metabolic pathway")
pathway_metabolite_only_list = df_pathway_metabolite_only["Pathway"].values.tolist()

df_pathways = df_pathways.drop("Features",axis=1)
df_pathways['Features'] = df_pathways['Pathways'].map(pathway_dict)
df_pathways['Genes'] = df_pathways['Pathways'].map(Gene_features)
df_pathways['Metabolites'] = df_pathways['Pathways'].map(Compound_features)
df_pathways['Reoccurrence'] = df_pathways['Pathways'].apply(lambda x: "Yes" if x in pathway_metabolite_only_list else "No")

# Save clean pathway data #
writer3 = pd.ExcelWriter(os.path.join(Folder2,File3_clean), engine='xlsxwriter')
df_pathways.to_excel(writer3,sheet_name='Jointed Pathway',index=False)
writer3.save()

## Jointed pathway results ##
# read results #
df_results = pd.read_csv(os.path.join(Folder1,File5))
df_results.rename(columns={"Unnamed: 0": "Pathways"},inplace=True)
df_results = df_results.replace({"Glycolysis or Gluconeogenesis":"Glycolysis / Gluconeogenesis"})
# save clean results #
writer4 = pd.ExcelWriter(os.path.join(Folder2,File4_clean), engine='xlsxwriter')
df_results.to_excel(writer4,sheet_name='Jointed Pathway results',index=False)
writer4.save()
