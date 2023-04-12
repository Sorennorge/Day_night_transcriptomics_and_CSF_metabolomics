# MacAulayLab Day-night fluctuations in choroid plexus transcriptomics and CSF metabolomics
The work and scripts are done by the MacAulay Lab.\
All programs used are free and open-source.
In the interest of open science and reproducibility, all data and source code used in our research is provided here.\
Feel free to copy and use code, but please cite:\
(coming soon) \
*Remember* rewrite file_names and folder_names suitable for your pipeline.\

## The RNAseq and metabolomics Analysis follows these steps:
## Raw data analysis - Library Build, Mapping and Quantification ##
The analysis uses RNA STAR for mapping and RSEM for TPM quantification.
### RNA-STAR and RSEM Library build and indexing ###

1.1.1 - RNA_STAR_Indexing.sh \
1.2.1 - RSEM_Indexing.sh

### RNA-STAR Mapping and RSEM quantification ###

1.1.2 -RNA_STAR_RNAseq2.sh \
1.2.2 - RSEM_RNAseq2.sh

## Metabolomics ##

### Data quantification and trimming ###

2.1.1 - Data analysis.py \

### PCA plots ###

2.1.2 - PCA plots.py \

### Pathway plot ###

2.2.1. - Pathway plot.py \

### Calculate padj (BH) ###

2.3.1 - Calculate adjusted pvalue.R \

### Volcano plot ###

2.4.1 - Volcano plot \

## RNA sequencing ##

### Raw counts table ###

3.1.1 - Raw data GeneCounts.py \
3.1.2 - CeneCount count table.py \
3.1.3 - Reduce GeneCount count table.py

### TPM count tables ###

3.1.4 - Raw data RSEM.py \
3.1.5 - RSEM Create count table.py \
3.1.6 - Reduce RSEM count table.py

### Differentially expressed genes ###

3.2.1 - DE Analysis Sleep.R

### PCA plots ###

3.2.2 - DE Analysis Sleep PCA.R

### Volcano plot ###

3.2.3 - DE Analysis Volcano.R

### Heatmap ###

3.2.4 - DE Analysis Heatmap.R

### Overview piechart ###

3.3.1 - Piechart of differentially expressed genes - Overview.py

### Enrichment analysis ###

3.4.1 - Enrichment analysis.py

### Create supplementary tables ###

4.1.1 - Supplementary tables.py

### Panther database analysis ###

5.1.1 - Panther DB data analysis.py

### Subcelluar location - Transport ###

5.1.2 - Subcellular location of transporters.py

### Transport analysis ###

5.1.3 - Transport analysis.py

## Jointed-pathway analysis ##

### Generate metaboanalyst input ###
website: https://metaboanalyst.ca/ \

6.1.0 - Generate jointed pathway input.py

### Clean data from metaboanalyst ###

6.1.1 - Clean Jointed pathway data.py

### Jointed-pathway analysis ###

6.1.2 - Jointed pathway analysis.py

### Jointed-pathway analysis - regulation ###

6.1.3 - Jointed pathway regulation analysis.py
