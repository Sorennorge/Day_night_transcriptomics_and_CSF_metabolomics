# -*- coding: utf-8 -*-

### Subcellular locations of transporters ###

## Libraries ##

import os
import pandas as pd

## Folders ##

Folder1 = "Data/RNA/Panther DB/Clean data"
Folder2 = "Data/RNA/Panther DB/Subcellular location"

## Files ##

File1 = "Transporter_list_correct.xlsx"
File2 = "Compartment_scores_of_Transporters.csv"
File3 = "Transporters_with_scores.xlsx"

## Load data ##

df_1 = pd.read_excel(os.path.join(Folder1,File1),sheet_name='Transporters')
df_2 = pd.read_csv(os.path.join(Folder2,File2),sep=',',na_values=("Yes"),engine="python")
df_2 = df_2.fillna(0.0)

df_3 = df_2[['query term','compartment::plasma membrane']]

df_3.rename(columns=({'query term':'Genes','compartment::plasma membrane':'Plasma membrane scores'}),inplace=True)

df_3.to_excel(os.path.join(Folder2,File3),sheet_name="Transport scores",index=False)
