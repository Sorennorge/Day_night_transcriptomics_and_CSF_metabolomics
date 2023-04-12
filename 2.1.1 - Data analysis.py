# -*- coding: utf-8 -*-
"""
Created on Mon Jan  2 10:28:47 2023

@author: dcs839
"""

### Metabolomics data analysis ###

## Libraries ##

import os
import numpy as np
import pandas as pd
import math
from scipy.stats import gmean
from scipy.stats import ttest_ind
from OUTLIERS import smirnov_grubbs as grubbs

## Functions ##

def log2FC_func(A,B):
    FC = B/A
    log2FC = round(math.log(FC,2),4)
    return log2FC

## Folders ##

Folder1 = "Data/Metabolomics/Raw data"
Folder2 = "Data/Metabolomics/Raw data analysis"
Folder3 = "Data/Metabolomics/Metabolomicsworkbench"
Folder4 = "Data/Metabolomics/PCA data"
Folder5 = "Data/Metabolomics/Data analysis"

os.makedirs(Folder2, exist_ok=True)
os.makedirs(Folder4, exist_ok=True)
os.makedirs(Folder5, exist_ok=True)

## Files ##

File1 = "Raw_data_transposed.csv"
File2 = "Raw_data_DP.xlsx"
File3 = "Metabolite_HMDB.xlsx"
File4 = "PCA_data.xlsx"
File5 = "PCA_targets.xlsx"
File6 = "Outlier overview.xlsx"
File7 = "DP overview.xlsx"
File8 = "Data analysis overview.xlsx"

## Global variables ##

# file load -> data type initialisation #

dtype_dict = {'Name':str}
for i in range(1,8,1):
    dtype_dict['QC{}'.format(i)] = float

for i in range(1,15,1):
    if i == 4 or i == 11:
        pass
    else:
        dtype_dict['N{}'.format(i)] = float
for i in range(1,12,1):
    dtype_dict['D{}'.format(i)] = float

# Dict variables #

Info_dict_day = {}
Info_dict_Night = {}
Info_dict_QC = {}
metabolite_list = []

## Raw data variables ##

# descriptive power #
Raw_DP_dict = {}

# Significant metabolite list (DP > 2.5)
Raw_Significant_metabolite_list = []

## Data analysis ##

# descriptive power calc #
DP_dict = {}
# Significant metabolite list (DP > 2.5)
Significant_metabolite_list = []

# Outliers #

Day_Outlier_dict = {}
Night_Outlier_dict = {}
Day_without_outliers = {}
Night_without_outliers = {}

# Normalized data dicts #

Normalized_Day_without_outliers = {}
Normalized_Night_without_outliers = {}

## mean, std, log2FC (mean), pvalue ##

Mean_day = {}
Mean_night = {}
Std_day = {}
Std_night = {}
log2fc_dict = {}
pvalue_dict = {}

### load data ###

df = pd.read_csv(os.path.join(Folder1,File1),delimiter=";",encoding = "cp1252",header=0,skiprows=[1],dtype=dtype_dict,decimal=",")

# Modulate data for dict conversion #

df_T = df.transpose()

data_dict = df_T.to_dict()

## Construct info dictionaries for calculations ##
for key in data_dict:
    # Metabolite list #
    metabolite_list.append(data_dict[key]['Name'])
    # Day #
    Day_array = []
    for i in range(1,12,1):
        Day_array.append(data_dict[key]['D{}'.format(i)])
    Info_dict_day[data_dict[key]['Name']] = np.array(Day_array,dtype=float)
    Night_array = []
    for i in range(1,15,1):
        if i == 4 or i == 11:
            pass
        else:
            Night_array.append(data_dict[key]['N{}'.format(i)])
    Info_dict_Night[data_dict[key]['Name']] = np.array(Night_array,dtype=float)
    QC_array = []
    for i in range(1,8,1):
        QC_array.append(data_dict[key]['QC{}'.format(i)])
    Info_dict_QC[data_dict[key]['Name']] = np.array(QC_array,dtype=float)
    
## Calculate DP for raw data and exclude metabolites for PCA plot ##

# calculate DP #
for key in metabolite_list:
    Raw_day_and_night_concatenated = np.concatenate((Info_dict_day[key], Info_dict_Night[key]), axis=None)
    Raw_DP = np.std(Raw_day_and_night_concatenated,ddof=1)/np.std(Info_dict_QC[key],ddof=1)
    Raw_DP_dict[key] = Raw_DP
    if Raw_DP > 2.5:
        Raw_Significant_metabolite_list.append(key)
    else:
        pass

# create dataframe from DP values (raw) #
raw_df_DP = pd.DataFrame({'Metabolites': list(Raw_DP_dict.keys()),'DP':list(Raw_DP_dict.values())})

raw_df_DP['Significance'] = np.where(raw_df_DP['DP'] >= 2.5, "Yes", "No")

if os.path.isfile(os.path.join(Folder2,File2)):
    pass
else:
    raw_df_DP.to_excel(os.path.join(Folder2,File2),sheet_name='Raw DP calc',header=True,index=False)

## Check if all Significant metabolites from raw data has HMDB numbers ##
# Running all metabolite names through metabolomcsworkbench.org to get valid HMDB numbers #
# Input file is the curated results from this workbench output #

HMDB_df = pd.read_excel(os.path.join(Folder3,File3))
HMDB_dict = HMDB_df.set_index('Metabolite').to_dict()['HMDB']
Raw_Significant_metabolite_list_HMDB = []

for key in HMDB_dict:
    if isinstance(HMDB_dict[key], str):
        if key in Raw_Significant_metabolite_list:
            Raw_Significant_metabolite_list_HMDB.append(key)
        else:
            pass
    else:
        pass

## Create PCA plot tables ##

# transpose data for correct PCA annotation #
PCA_df_T = df.set_index('Name').T
# Only include metabolites with DP > 2.5 #
PCA_df = PCA_df_T.loc[:, PCA_df_T.columns.isin(Raw_Significant_metabolite_list_HMDB)]

PCA_targets = pd.DataFrame(index=PCA_df.index.copy())
PCA_targets['Samples'] = PCA_targets.index
PCA_targets.loc[PCA_targets['Samples'].str.contains('QC'), 'Target'] = 'QC'
PCA_targets.loc[PCA_targets['Samples'].str.contains('D'), 'Target'] = 'Day'
PCA_targets.loc[PCA_targets['Samples'].str.contains('N'), 'Target'] = 'Night'

if os.path.isfile(os.path.join(Folder4,File4)):
    pass
else:
    PCA_df.to_excel(os.path.join(Folder4,File4),sheet_name='PCA data',header=False,index=False)

if os.path.isfile(os.path.join(Folder4,File5)):
    pass
else:
    PCA_targets.to_excel(os.path.join(Folder4,File5),sheet_name='PCA Targets',header=False,index=False)

## For PCA plot ##
# Run script 1.2.1 - PCA plot.py #

### Data analysis ###

## Outlier testting ##
for key in metabolite_list:
    Day_Outlier_dict[key] = []
    Night_Outlier_dict[key] = []
    Day_without_outliers[key] = grubbs.test(Info_dict_day[key], alpha=.05)
    Night_without_outliers[key] = grubbs.test(Info_dict_Night[key], alpha=.05)
    
    outlier_index_Day = grubbs.two_sided_test_indices(Info_dict_day[key], alpha=.05)
    outlier_index_Night = grubbs.two_sided_test_indices(Info_dict_Night[key], alpha=.05)
    for i in outlier_index_Day:
        outlier_sample_day = "D{}".format(i+1)
        Day_Outlier_dict[key].append(outlier_sample_day)
    for i in outlier_index_Night:
        if i <= 2:
            outlier_sample_Night = "N{}".format(i+1)
        elif i > 2 and i < 9:
            outlier_sample_Night = "N{}".format(i+2)
        else:
            outlier_sample_Night = "N{}".format(i+3)
        Night_Outlier_dict[key].append(outlier_sample_Night)

# create outlier dataframe and save outlier file #

df_outlier = pd.DataFrame({'Metabolites': list(Day_Outlier_dict.keys()),'Outliers day':list(Day_Outlier_dict.values()),'Outliers night':list(Night_Outlier_dict.values())})
df_outlier['Outliers day'] = df_outlier['Outliers day'].str.join(",")
df_outlier['Outliers night'] = df_outlier['Outliers night'].str.join(",")

df_outlier.to_excel(os.path.join(Folder5,File6),sheet_name='Outliers',header=True,index=False)


## calculate DP and assign significant metabolites after outliers have been removed ##
for key in metabolite_list:
    Day_and_night_concatenated = np.concatenate((Day_without_outliers[key], Night_without_outliers[key]), axis=None)
    DP = np.std(Day_and_night_concatenated,ddof=1)/np.std(Info_dict_QC[key],ddof=1)
    DP_dict[key] = DP
    if DP > 2.5:
        Significant_metabolite_list.append(key)
    else:
        pass
# Create DP dataframe and save metabolite and descriptive power (DP) to file #
df_DP = pd.DataFrame({'Metabolites': list(DP_dict.keys()),'DP':list(DP_dict.values())})
df_DP['Significance'] = np.where(df_DP['DP'] >= 2.5, "Yes", "No")
df_DP.to_excel(os.path.join(Folder5,File7),sheet_name='Descriptive power',header=True,index=False)
 
## excluded significant metabolites without a HMDB number ##
Significant_metabolites_list_HMDB = [] 
for key in HMDB_dict:
    if isinstance(HMDB_dict[key], str):
        if key in Significant_metabolite_list:
            Significant_metabolites_list_HMDB.append(key)
        else:
            pass
    else:
        pass

## Normalize data based on geomean of the QC samples ##

for key in Significant_metabolites_list_HMDB:
    geometric_mean_QC = gmean(Info_dict_QC[key])
    Normalized_Day_without_outliers[key] = Day_without_outliers[key]/geometric_mean_QC
    Normalized_Night_without_outliers[key] = Night_without_outliers[key]/geometric_mean_QC

## calculate mean, std, log2FC (mean), pvalue ##

for key in Significant_metabolites_list_HMDB:
    # mean #
    mean_day = np.mean(Normalized_Day_without_outliers[key])
    mean_night = np.mean(Normalized_Night_without_outliers[key])
    Mean_day[key] = mean_day
    Mean_night[key] = mean_night
    # Standard diviation #
    std_day = np.std(Normalized_Day_without_outliers[key],ddof=1)
    std_night = np.std(Normalized_Night_without_outliers[key],ddof=1)
    Std_day[key] = std_day
    Std_night[key] = std_night
    # Log2FC for mean #
    log2fc_mean = log2FC_func(mean_day,mean_night)
    log2fc_dict[key] = log2fc_mean
    # Welch t-test without equal variance
    welch_p_value = ttest_ind(Normalized_Day_without_outliers[key], Normalized_Night_without_outliers[key], equal_var = False)[1]
    pvalue_dict[key] = welch_p_value
## Create overview dataframe ##    
df_overview = pd.DataFrame({'Metabolites': list(Significant_metabolites_list_HMDB)})
# add HMDB #
df_overview['HMDB'] = df_overview['Metabolites'].map(HMDB_dict) 
# Add mean and std #
df_overview['Day (mean)'] = df_overview['Metabolites'].map(Mean_day)
df_overview['Day (std)'] = df_overview['Metabolites'].map(Std_day)
df_overview['Night (mean)'] = df_overview['Metabolites'].map(Mean_night)
df_overview['Night (std)'] = df_overview['Metabolites'].map(Std_night)
# Add Log2FC and p-value #
df_overview['Log2FC'] = df_overview['Metabolites'].map(log2fc_dict)
df_overview['P-value'] = df_overview['Metabolites'].map(pvalue_dict)

# Save overview to excel file #
df_overview.to_excel(os.path.join(Folder5,File8),sheet_name='Overview',header=True,index=False)
