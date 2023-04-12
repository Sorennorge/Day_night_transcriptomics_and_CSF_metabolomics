# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 08:48:38 2023

@author: dcs839
"""

### Piechart - differentially expressed overview ###

## Libraries ##

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

    
## Files ##

File1 = "Count_table_Sleep_RSEM_Reduced.csv"

File1_out = "Differentially expressed piechart plot.png"
File2_out = "Differentially expressed piechart plot (clean).png"

## Folders ##

Folder1 = "Data/RNA/Count Tables"
Folder2 = "Results/RNA/Piecharts"

os.makedirs(Folder2, exist_ok=True)

## Variables ##

# From Deseq2 #

up_regulated = 1496
down_regulated = 1282
Count_deseq2 = 21025-6473-143 # all - lowcount - outliers

# Percentage
Deseq_up = round(up_regulated/Count_deseq2*100,1)
Deseq_down = round(down_regulated/Count_deseq2*100,1)
Deseq_non = round(100-Deseq_up-Deseq_down,1)

# Set label data
Labels_1 = ['Not differentially expressed ({} %)'.format(Deseq_non),
            'Downregulated ({} %)'.format(Deseq_down),
            'Upregulated ({} %)'.format(Deseq_up)]
# Set data
data_1 = [Deseq_non,Deseq_down,Deseq_up]

# Set colors #
colors_set = ['#B3B3B3','#00bfff','#ee2c2c']

# Plot piechart of differentially expressed genes #
plt.figure(figsize=(20,20))
plt.pie(data_1, labels = Labels_1, colors = colors_set,explode=[0,0.2,0.2],textprops={'fontsize': 60},wedgeprops = {"edgecolor":"black",'linewidth': 3},startangle=60,pctdistance=0.8,labeldistance=1.1)
plt.title('Differentially expressed genes - overview ', fontsize=60)
plt.savefig(os.path.join(Folder2,File1_out),dpi=600,bbox_inches='tight')
#plt.show()

#Plot piechart without labels
plt.figure(figsize=(20,20))
plt.pie(data_1, colors = colors_set,explode=[0,0.2,0.2],textprops={'fontsize': 60},wedgeprops = {"edgecolor":"black",'linewidth': 3},startangle=60,pctdistance=0.8,labeldistance=1.1)
plt.savefig(os.path.join(Folder2,File2_out),dpi=600,bbox_inches='tight',transparent=True)
#plt.show()

