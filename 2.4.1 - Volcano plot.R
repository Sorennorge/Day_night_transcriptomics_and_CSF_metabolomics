### Metabolomics - Volcano plot ###
# The packages can be installed within R with:
# intstall.packages(c("tidyverse", "ggrepel"))

## Libraries ##
library(tidyverse)
library("readxl")
library("openxlsx")

# set working directory #
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Folders #

Folder1 <- "Results/Metabolomics/Tables"
Folder2 <- "Results/Metabolomics/Volcano"

dir.create(file.path(Folder2), showWarnings = FALSE)

# Files #

File1 <- "Data analysis overview version 2.xlsx"
File2 <- "Volcano plot.png"

# Read data #

data <- read_excel(file.path(Folder1, File1),sheet = "Overview") 

# Change Adjusted P-value to padj for data manipulation #
colnames(data)[8] <- "Pvalue"
colnames(data)[9] <- "Padj"

# mutate data (annotate up/down regulation)

data <- data %>% 
  mutate(
    Expression = case_when(Log2FC > 0.0 & Padj <= 0.1 ~ "Up-regulated",
                           Log2FC < 0.0 & Padj <= 0.1 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )


# Create volcano plot #

p2 <- ggplot(data, aes(Log2FC , -log(Padj,10))) +
  geom_hline(yintercept=-log10(0.1),linetype = "dashed", col="grey50")+
  geom_point(aes(fill = Expression), size = 2,shape=21)+
  scale_x_continuous(name=expression("log"[2]*"FC"),limits = c(-2,2))+
  scale_y_continuous(name=expression("-log"[10]*"pvalue"),limits = c(0,1.7))+
  scale_fill_manual(values = c("dodgerblue3", "gray80", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  theme_bw()+
  theme(legend.position="none")
# plot #
p2

# Create and save volcano plot as file #
png(filename = file.path(Folder2, File2),
    width = 12, height = 8, units = "in",
    bg = "transparent",res=800)
p2
dev.off()
