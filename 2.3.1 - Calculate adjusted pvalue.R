### Metabolomics ###

## Calculate adjusted p-value ##

# Libraries #
library("readxl")
library("openxlsx")

# set working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Folders #

Folder1 <- "Data/Metabolomics/Data analysis"
Folder2 <- "Results/Metabolomics/Tables"

dir.create(file.path(Folder2), showWarnings = FALSE)

# Files #

File1 <- "Data analysis overview.xlsx"
File2 <- "Data analysis overview version 2.xlsx"

# Read data #

data <- read_excel(file.path(Folder1, File1),sheet = "Overview") 

# Calculate adjusted p-value (padj) by Benjamini and Hochberg method #

data$'Adjusted P-value' <- p.adjust(data$'P-value',method="BH")

# Save table to file #
write.xlsx(data,file.path(Folder2, File2),colNames = TRUE,sheetName="Overview")

