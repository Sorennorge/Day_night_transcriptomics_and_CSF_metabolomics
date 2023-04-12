### If DEseq2 needs to be insalled ###
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.14")

library("DESeq2")
library("ggplot2")
library("tidyverse")
library("cowplot")

## set working directory ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Folders

Folder1 = "Data/RNA/Count Tables"
Folder2 = "Results/RNA/Volcano"

# check if sub directory exists 
if (file.exists(Folder2)){
  invisible()
} else {
  # create a directory if it doesnt exists
  dir.create(Folder2,recursive = TRUE)
}

## Files ##

File1 = "Count_table_Sleep_RNASTAR_Reduced.csv"
File2 = "Volcano plot.png"

### load data ###
data <- read.table(file.path(Folder1,File1), header=T, sep=";", stringsAsFactors=F, row.names=1)

## create condition / treatment ##
sample_info <- data.frame(condition = factor(c(rep("Day",6),c(rep("Night", 6)))),row.names = factor(colnames(data)))

## create DE dataset in DESeq2 ##

DESeq.ds <- DESeqDataSetFromMatrix(countData = round(data), colData = sample_info, design = ~ condition)
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0 , ]
colData(DESeq.ds)$treatment <- relevel(colData(DESeq.ds)$condition, "Night")

### Differential expression analysis with DESeq2 ##

DESeq.ds <- DESeq(DESeq.ds)
DGE.results <- results(DESeq.ds, independentFiltering = TRUE, alpha = 0.05)
DGE.results.sorted <- DGE.results[order(DGE.results$log2FoldChange), ]

## Assign differentially expressed labels ##

# Unchanged -> not differentially expressed #
DGE.results.sorted$diffexpressed <- "Unchanged"
# differentially expressed (up) #
DGE.results.sorted$diffexpressed[DGE.results.sorted$log2FoldChange > 0.0 & DGE.results.sorted$padj < 0.05] <- "Up-regulated"
# differentially expressed (down) #
DGE.results.sorted$diffexpressed[DGE.results.sorted$log2FoldChange < 0.0 & DGE.results.sorted$padj < 0.05] <- "Down-regulated"

## Sort the data so the differential expressed is on top ##
data_volcano <- data.frame(x = DGE.results.sorted$log2FoldChange, y = DGE.results.sorted$padj, diff = DGE.results.sorted$diffexpressed)
data_volcano.sorted = data_volcano[order(data_volcano$diff), ]
# Remove NA rows from dataset #
data_volcano.sorted <- data_volcano.sorted[complete.cases(data_volcano.sorted$y),]

## Choose color scheme ##
mycolors <- c("deepskyblue1", "firebrick2", "grey70")
names(mycolors) <- c("Down-regulated", "Up-regulated", "Unchanged")

# If adjusted pvalues are extremely -> set to max 10e-250
data_volcano.sorted <- data_volcano.sorted %>% 
  mutate(
    y2 = case_when(y <= 1e-200 ~ 1e-200,
                                 TRUE ~ y)
  )

## GGplot ##
p <- ggplot(data_volcano.sorted, aes(x,-log10(y2),colour=diff)) + 
  geom_point(aes(colour=diff, fill=diff),size=1,shape=21,colour = "black",stroke = 0.2) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = mycolors) +
  scale_x_continuous(name=expression("log"[2]*"FC"), limits=c(-5.1, 4.1),breaks=seq(-5,4,1)) +
  scale_y_continuous(name=expression("-log"[10]*"pvalue"))+
  theme(legend.position="none")

p

# In our case, the values are extremely high, so we'll be adjusting the layout to fix the mix #


## Creating two graphs with different y axis ##

## Plot volcano ##
p2 <- ggplot(data_volcano.sorted, aes(x,-log10(y2),colour=diff)) + 
  geom_point(aes(colour=diff, fill=diff),size=2,shape=21,colour = "black",stroke = 0.2) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = mycolors) +
  scale_x_continuous(name=expression("log"[2]*"FC"), limits=c(-5.1, 4.1),breaks=seq(-5,4,1)) +
  scale_y_continuous(name=expression("-log"[10]*"(adjusted pvalue)"),limits = c(0,20),breaks = seq(0,20, 5))+
  theme(legend.position="none")
p2


p3 <- ggplot(data_volcano.sorted, aes(x,-log10(y2),colour=diff)) + 
  geom_point(aes(colour=diff, fill=diff),size=2,shape=21,colour = "black",stroke = 0.2) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_fill_manual(values = mycolors) +
  scale_x_continuous(name=NULL, limits=c(-5.1, 4.1),breaks=seq(-5,4,1),labels=NULL) +
  scale_y_continuous(name=NULL,limits = c(20,200),breaks = c(20,50,100,150,200))+
  theme(legend.position="none")
p3

## Using cowplot to add these together ##
  
p4 <- plot_grid(p3, p2, align = "v", nrow = 2,rel_heights = c(3/10, 7/10))
p4

## Save volcano plot as image file ##
ggsave(file.path(Folder2,File2),plot=p4,width=6,height=5,units=c("in"),dpi=800)
