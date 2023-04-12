# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:33:16 2023

@author: dcs839
"""

### Jointed pathway regulation ###

## Libs ##

import os
import pandas as pd
import numpy as np
import math

## Folders ##

Folder1 = "Data/Metabolomics/Jointed pathway analysis data/Clean data"
Folder2 = "Data/Metabolomics/Jointed pathway analysis data/Metabolyst input"
Folder3 = "Results/Jointed pathway"

## Files ##

File1 = "Metabolite_name_mapping.xlsx"
File2 = "Gene_name_mapping.xlsx"
File3 = "Jointed_Pathway_data.xlsx"

File4 = "Metabolomics_compound_logFC.txt"
File5 = "RNA_Gene_symbol_logFC.txt"

File6 = "Regulation overview.xlsx"

## Load data ##

df_1 = pd.read_excel(os.path.join(Folder1,File1),dtype=str)
df_2 = pd.read_excel(os.path.join(Folder1,File2),dtype=str)
df_3 = pd.read_excel(os.path.join(Folder1,File3),dtype=str)
df_3 = df_3[df_3['Reoccurrence'] == 'Yes']

df_4 = pd.read_csv(os.path.join(Folder2,File4),sep="\t")
df_5 = pd.read_csv(os.path.join(Folder2,File5),sep="\t")

## Variables ##

from_Compound_to_HMDB = df_1.set_index('Compound')['HMDB'].to_dict()
from_Entrez_to_Gene = df_2.set_index('Entrez')['Gene'].to_dict()

Pathway_genes = df_3.set_index('Pathways')['Genes'].to_dict()
Pathway_compounds = df_3.set_index('Pathways')['Metabolites'].to_dict()

LogFC_compound = df_4.set_index('#compound')['log2FC'].to_dict()
LogFC_genes = df_5.set_index('#Official')['log2FC'].to_dict()

## set pathway genes to lists ##

Pathway_genes_remapped = {}
for key in Pathway_genes:
    if isinstance(Pathway_genes[key], str):
        gene_list = Pathway_genes[key].split("|")
        gene_list_remapped = [from_Entrez_to_Gene[x] for x in gene_list]
        Pathway_genes_remapped[key] = gene_list_remapped
    else:
        Pathway_genes_remapped[key] = []

Pathway_genes_remapped_FC = {}
for key in Pathway_genes_remapped:
    Pathway_genes_remapped_FC[key] = np.array([LogFC_genes[x] for x in Pathway_genes_remapped[key]])

Pathway_gene_stats = {}
Pathway_gene_stats_total = {}
Pathway_gene_stats_up = {}
Pathway_gene_stats_down = {}
for key in Pathway_genes_remapped_FC:
    Number_total = np.size(Pathway_genes_remapped_FC[key])
    Upreg = np.count_nonzero(Pathway_genes_remapped_FC[key] > 0)
    Downreg = np.count_nonzero(Pathway_genes_remapped_FC[key] < 0)
    Pathway_gene_stats[key] = [Number_total,Upreg,Downreg]
    Pathway_gene_stats_total[key] = Number_total
    Pathway_gene_stats_up[key] = Upreg
    Pathway_gene_stats_down[key] = Downreg

Pathway_compounds_remapped = {}
for key in Pathway_compounds:
    if isinstance(Pathway_compounds[key], str):
        compound_list = Pathway_compounds[key].split("|")
        compound_list_remapped = [from_Compound_to_HMDB[x] for x in compound_list]
        Pathway_compounds_remapped[key] = compound_list_remapped
    else:
        Pathway_compounds_remapped[key] = []    

Pathway_compounds_remapped_FC = {}
for key in Pathway_compounds_remapped:
    Pathway_compounds_remapped_FC[key] = np.array([LogFC_compound[x] for x in Pathway_compounds_remapped[key]])

Pathway_compound_stats_total = {}
Pathway_compound_stats_up = {}
Pathway_compound_stats_down = {}
for key in Pathway_compounds_remapped_FC:
    Number_total = np.size(Pathway_compounds_remapped_FC[key])
    Upreg = np.count_nonzero(Pathway_compounds_remapped_FC[key] > 0)
    Downreg = np.count_nonzero(Pathway_compounds_remapped_FC[key] < 0)
    Pathway_compound_stats_total[key] = Number_total
    Pathway_compound_stats_up[key] = Upreg
    Pathway_compound_stats_down[key] = Downreg

df_overview = df_3['Pathways']
df_overview = df_overview.to_frame()
df_overview['Compounds (HMDB)'] = df_overview['Pathways'].map(Pathway_compounds_remapped)
df_overview['Compounds (log2FC)'] = df_overview['Pathways'].map(Pathway_compounds_remapped_FC)
df_overview['Genes'] = df_overview['Pathways'].map(Pathway_genes_remapped)
df_overview['Genes (log2FC)'] = df_overview['Pathways'].map(Pathway_genes_remapped_FC)
df_overview['Compounds (Total)'] = df_overview['Pathways'].map(Pathway_compound_stats_total)
df_overview['Compounds (Upreg)'] = df_overview['Pathways'].map(Pathway_compound_stats_up)
df_overview['Compounds (Downreg)'] = df_overview['Pathways'].map(Pathway_compound_stats_down)
df_overview['Genes (Total)'] = df_overview['Pathways'].map(Pathway_gene_stats_total)
df_overview['Genes (Upreg)'] = df_overview['Pathways'].map(Pathway_gene_stats_up)
df_overview['Genes (Downreg)'] = df_overview['Pathways'].map(Pathway_gene_stats_down)

writer = pd.ExcelWriter(os.path.join(Folder3,File6), engine='xlsxwriter')
df_overview.to_excel(writer,sheet_name='Regulation overview',index=False)
writer.save()