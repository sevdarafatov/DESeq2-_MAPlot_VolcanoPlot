library(magrittr)
library(dplyr)
library(readr)
library(DESeq2)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(tibble)
library(pheatmap)
library(tidyr)
library(reshape2)
library(VennDiagram)
library(gprofiler2)
library(debrowser)
library('biomaRt')
library("org.Hs.eg.db")
library(tidyverse) # already know about this from Step 1 script
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot)
library("ggrepel")
# load the file


df2 <- read.csv("feature_counts9_result.csv", header=TRUE, row.names = 1, sep = '\t' )
head(df2)

#filter counts and remove less than 10
counts <- df2 %>%  filter_all(all_vars(. > 10))
condition <- c( "control", "control","control","case", "case","case")
metadata <- data.frame(condition)
rownames(metadata) <- colnames(counts[2:7])

#DEseq2
dds_smoc2 <- DESeqDataSetFromMatrix(countData = counts[2:7],
                                    colData = metadata,
                                    design = ~ condition)
heat_colors <- brewer.pal(n = 6, name = "Paired")

# Determine the size factors to use for normalization
dds_smoc2 <- estimateSizeFactors(dds_smoc2)

# Extract the counts
smoc2_counts <- counts(dds_smoc2)

# Extract the normalized counts
smoc2_normalized_counts <- counts(dds_smoc2, normalized=TRUE)

#boxplot
boxplot(log2(smoc2_counts), ylim = c(0, 20), outline=TRUE,main="Log2 Counts per Million (CPM)" ,sub="filtered, non-normalised" ,ylab="Log2 expression", xlab="Samples", col= heat_colors, outer=F)
boxplot(log2(smoc2_normalized_counts), ylim = c(0, 20), outline=TRUE,main="Log2 Counts per Million (CPM)", sub="filtered, normalised", ylab="Log2 expression",xlab="Samples", col= heat_colors)

# density plot (representation of the distribution of a numeric variable)
LOG <- data.frame(log2(smoc2_normalized_counts))
da <- tibble::rowid_to_column(LOG, "ID")
molted <- melt(da,id.vars="ID")
molted%>%
  ggplot(aes(x=value, fill=variable)) +
  geom_density(alpha=0.3) +
  labs(x= "Log2 expression", y= "Density") +
  scale_fill_manual('Samples',values=heat_colors)

# Transform the normalized counts
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Plot the PCA of PC1 and PC2
plotPCA(vsd_smoc2)
dds_smoc2 <- DESeq(dds_smoc2)

# Extract the results of the differential expression analysis between two groups only
smoc2_res <- results(dds_smoc2,
                     contrast = c("condition", "case", "control"),
                     alpha = 0.05)
#clean data
smoc2_res_all <- data.frame(smoc2_res)
all <- data.frame(smoc2_res_all) %>% mutate(threshold = padj < 0.05)
all$sig <- ifelse(all$threshold == 'FALSE', 'FALSE', ifelse(all$threshold == 'TRUE' & all$stat > 0 , 'Up', 'Down'))
all$threshold <- c()
all$Genes <-ifelse(rownames(all) == rownames(counts) ,counts$GeneSymbol,'')
res <- all %>% drop_na()

unique(res$Genes)
# Create MA plot and main has the groups
up <- arrange(subset(res, sig == "Up"), desc(log2FoldChange))
down <- arrange(subset(res, sig == "Down"), desc(log2FoldChange))
nonsig <- arrange(subset(res, sig == "FALSE"), desc(log2FoldChange))
with(nonsig, plot(y=log2FoldChange, baseMean, pch = 3,  main="MAplot",col="grey", cex = 0.2, xlim=c(0,10000), ylim=c(-6,6)))+
  abline(h = 0, col = "grey")
with(subset(up, sig =="Up"), points(y=log2FoldChange,baseMean, pch = 3, cex = 0.2, col = "blue"))
with(subset(down, sig=='Down'), points(y=log2FoldChange,baseMean, pch = 3,cex = 0.2, col = "red"))
legend("bottomright", inset=.01,c("Up", "Down" , "Non-SIG"), fill= c("blue","red","grey"), horiz=F, cex=0.8)
selected <- c(1:5)
text(log2FoldChange[selected] ~ baseMean[selected], labels= Genes[selected], data=up,
     pos = 3, offset =0.4, cex = 0.5, col = "black", font = 1)
selected <- c(nrow(down), nrow(down)-1, nrow(down)-2, nrow(down)-3,nrow(down)-4)
text(log2FoldChange[selected] ~ baseMean[selected], labels= Genes[selected], data=down,
     pos = 1, offset =0.4, cex = 0.5, col = "black", font = 1)

# Create the Volcano plot
res$log10padj <- -log10(res$padj)
upv <- arrange(subset(res, sig == "Up"), desc(log10padj))
downv <- arrange(subset(res, sig == "Down"), desc(log10padj))
nonsigv <- arrange(subset(res, sig == "FALSE"), desc(log10padj))
with(nonsigv, plot(y=log10padj, x= log2FoldChange, pch = 3,  main="Volcano plot",col= 'grey', cex = 0.2, xlim=c(-10,10), ylim=c(0,20)))
with(subset(upv), points(y=log10padj, x= log2FoldChange, pch = 3,cex = 0.2, col = "blue"))
with(subset(downv), points(y=log10padj, x= log2FoldChange, pch = 3,cex = 0.2, col = "red"))
legend("bottomright", inset=.01,c("Up", "Down" , "Non-SIG"), fill= c("blue","red","grey"), horiz=F, cex=0.8)
selected <- c(1,2,3,4,5)
text(log10padj[selected] ~ log2FoldChange[selected], labels= Genes[selected], data=upv,
     pos = 4, offset =0.2, cex = 0.5, col = "black", font = 1)
text(log10padj[selected] ~ log2FoldChange[selected], labels= Genes[selected], data=downv,
     pos = 4, offset =0.2, cex = 0.5, col = "black", font = 1)
     
     
