### If DEseq2 needs to be insalled ###
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.14")

library("DESeq2")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")

## set working directory ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Folders

Folder1 = "Data/RNA/Count Tables"
Folder2 = "Results/RNA/Heatmap"

# check if sub directory exists 
if (file.exists(Folder2)){
  invisible()
} else {
  # create a directory if it doesnt exists
  dir.create(Folder2,recursive = TRUE)
}


## Files ##

File1 = "Count_table_Sleep_RNASTAR_Reduced.csv"
File2 = "Heatmap.png"
File3 = "Heatmap_matrix_zscore.csv"

### load data ###
data <- read.table(file.path(Folder1,File1), header=T, sep=";", stringsAsFactors=F, row.names=1)

## create condition / treatment ##
sample_info <- data.frame(condition = factor(c(rep("Day",6),c(rep("Night", 6)))),row.names = factor(colnames(data)))

## create DE dataset in DESeq2 ##

DESeq.ds <- DESeqDataSetFromMatrix(countData = round(data), colData = sample_info, design = ~ condition)
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0 , ]
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, "Day")
DESeq.ds <- estimateSizeFactors(DESeq.ds)
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE)

## Create rlog and log norm counts ##
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE )
rlog.norm.counts <- assay(DESeq.rlog)
log.norm.counts <- log2(counts.sf_normalized +1)

### Differential expression analysis with DESeq2 ##
colData(DESeq.ds)$treatment <- relevel(colData(DESeq.ds)$condition, "Night")
DESeq.ds <- DESeq(DESeq.ds)
DGE.results <- results(DESeq.ds, independentFiltering = TRUE, alpha = 0.05)
summary(DGE.results)

DGE.results.sorted <- DGE.results[order(-DGE.results$log2FoldChange), ]
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))

heatmat_DGEgenes <- log.norm.counts[DGEgenes,]

### Set a color palette
heat_colors <- brewer.pal(10, "RdYlBu")

### Run pheatmap

p<-pheatmap(heatmat_DGEgenes, 
         color = rev(heat_colors),
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         border_color = "black", 
         fontsize = 5, 
         scale = "row", 
         fontsize_row = 5, 
         height = 20)
p

## Save plot to file ##

ggsave(file.path(Folder2,File2),plot=p,width=8,height=6,units=c("in"),dpi=800)

## Get the Zscores of the matrix ##

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
zscores <- scale_rows(heatmat_DGEgenes)

write.csv2(zscores,file.path(Folder2, File3))
