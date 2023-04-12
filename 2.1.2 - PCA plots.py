# -*- coding: utf-8 -*-

import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set()

## Folders ##

Folder1 = "Data/Metabolomics/PCA data"
Folder2 = "Results/Metabolomics/PCA"

os.makedirs(Folder2, exist_ok=True)

## Files ##

file_data = "PCA_data.xlsx"
file_targets = "PCA_targets.xlsx"

File_out = "PCA_plot.png"

## Load data ##
df_data = pd.read_excel(os.path.join(Folder1,file_data),header=None)

df_targets = pd.read_excel(os.path.join(Folder1,file_targets),header=None)

# data scaling
x_scaled = StandardScaler().fit_transform(df_data)

# set principal compnents #
n = 10
pca = PCA(n_components=n)

# transport data
pca_features = pca.fit_transform(x_scaled)

columns_n = []
for i in range(1,n+1,1):
    columns_n.append("PC{}".format(i))
# create dataframe with the n PC
pca_df = pd.DataFrame(
    data=pca_features, 
    columns=columns_n)

## Variance for plot 

Variance_PC_array = pca.explained_variance_ratio_

PC_1_V = round(Variance_PC_array[0]*100,2)
PC_2_V = round(Variance_PC_array[1]*100,2)

# map target names to PCA features   

pca_df['Target'] = df_targets[1]


plt.figure(figsize=(20,20))
sns.lmplot(
    x='PC2', 
    y='PC1', 
    data=pca_df, 
    hue='Target', 
    fit_reg=False, 
    legend=True
    )

plt.title('PCA Day and Night')
plt.xlabel( "PC 2 ({} %)".format(PC_2_V))
plt.ylabel( "PC 1 ({} %)".format(PC_1_V))
plt.ylim(-21,21)
plt.xlim(-10,10)
plt.savefig(os.path.join(Folder2,File_out),dpi=800,bbox_inches='tight')
plt.show()
