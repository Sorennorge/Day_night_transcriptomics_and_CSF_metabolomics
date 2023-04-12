### If DEseq2 needs to be insalled ###
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.14")

## Import libraries ##
library("DESeq2")
library("ggplot2")
library("openxlsx")

## set working directory ##
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Folders

Folder1 = "Data/RNA/Count Tables"
Folder2 = "Results/RNA/PCA"

# check if sub directory exists 
if (file.exists(Folder2)){
  invisible()
} else {
  # create a directory if it doesnt exists
  dir.create(Folder2,recursive = TRUE)
}


## Files ##

File1 = "Count_table_Sleep_RNASTAR_Reduced.csv"
File2 = "PCA data.xlsx"
File3 = "PCA plot.png"

## Dataset ##
data <- read.table(file.path(Folder1,File1), header=T, sep=";", stringsAsFactors=F, row.names=1)

# create condition / treatment #
sample_info <- data.frame(condition = factor(c(rep("Day",6),c(rep("Night", 6)))),row.names = factor(colnames(data)))

## create DE dataset in DESeq2 ##

DESeq.ds <- DESeqDataSetFromMatrix(countData = round(data), colData = sample_info, design = ~ condition)
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0 , ]
colData(DESeq.ds)$treatment <- relevel(colData(DESeq.ds)$condition, "Night")

## Create rlog and log norm counts ##


## Differential expression analysis with DESeq2 ##

DESeq.ds <- DESeq(DESeq.ds)
DGE.results <- results(DESeq.ds, independentFiltering = TRUE, alpha = 0.05)
summary(DGE.results)

## Create rlog dataframe ##
rlt <- DESeq2::rlogTransformation(DESeq.ds)

## Create PCA plot ##
p1 <- plotPCA(rlt, intgroup="condition", returnData=TRUE)

## Plot PCA ##
percentVar <- round(100 * attr(p1, "percentVar"))
p2 <- ggplot(p1,aes(PC1, PC2))+
  geom_point(aes(color=condition, fill=condition),size=5,color="black",shape=21,stroke = 1)+
  scale_fill_manual(values=c("#fdd835","#515151"))+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),legend.title = element_blank())
p2

## Write PCA xy coordinates to file ##

write.xlsx(p1,file.path(Folder2, File2),colNames = TRUE,sheetName="PCA")

## Save PCA plot as image ##

ggsave(file.path(Folder2, File3),plot=p2,width=8,height=6,units=c("in"),dpi=800)
