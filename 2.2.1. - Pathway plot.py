# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 14:32:21 2022

@author: dcs839
"""

### Pathway plot ###

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="white")
# folders #

Folder1 = "Data/Metabolomics/Pathway plot data"
Folder2 = "Results/Metabolomics/Fig1 Pathway"

os.makedirs(Folder2, exist_ok=True)

# files #

file_in = "Overview.xlsx"
file_out = "Metabolic Pathway.png"


df = pd.read_excel(os.path.join(Folder1,file_in),sheet_name="Pathway stats")
df2 = df[df['FDR'] < 0.1]
n = 20
diverging_colors = sns.color_palette("mako_r", n)

with sns.plotting_context(rc={"legend.fontsize":14}):
    plt.figure(figsize=(25,15))
    plt.tight_layout()
    ax = sns.relplot(data=df2, x="Impact", y="Ratio metab", hue="Group",palette=[diverging_colors[2],diverging_colors[int(n/2)]],size="#",sizes=(50, 400),edgecolor="black",linewidth=0.75)
    plt.title("Metabolic Pathway", size=20)
    h,l = ax.get_legend_handles_labels()
    plt.legend(h[0:3],l[0:3],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=13)
    plt.savefig(os.path.join(Folder2,file_out),dpi=800,bbox_inches='tight')
    plt.show()
