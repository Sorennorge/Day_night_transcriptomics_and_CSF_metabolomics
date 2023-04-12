# -*- coding: utf-8 -*-

### Enrichment plots ###

## Libraries ##

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## Folders ##

Folder1 = "Data/RNA/Enrichment data"
Folder2 = "Results/RNA/Enrichement Analysis"

os.makedirs(Folder2, exist_ok=True)
    
## Files ##

# The enrichment data are either constructed by using the panther database and manually constructing the data #
# Or use script 4.1.1 Panther db data analysis #
File1 = "Enrichement_all_DE_protein_classes_reduced.csv"
File2 = "Enrichment plot.png"

df_all = pd.read_csv(os.path.join(Folder1,File1),sep=";")

data_all = df_all['Value']
labels_all = df_all['Group']

n=60
colors = sns.color_palette('blend:#FFFFFF,#777696,#000000',n_colors=n)
colors[22] = colors[10]
color_set = ["#9696b4"]*15
color_set[2] = "#777696"
explode_set = [0.0]*15
explode_set[2] = 0.15

plt.figure(figsize=(20,20))
plt.pie(data_all, labels = labels_all, colors = colors[20:],explode=explode_set, autopct='%.0f%%',textprops={'fontsize': 30},wedgeprops = {"edgecolor":"black",'linewidth': 1.2},counterclock=False,startangle=90,pctdistance=0.8,labeldistance=1.1)
plt.savefig(os.path.join(Folder2,File2),dpi=600,bbox_inches='tight')
#plt.show()
