# BRCA1, BRCA2, WT groups
library('DESeq2') # DESeq2_1.26.0
directory <- # Directory in wich Kallisto tsv files exist
fileName<- c("T_BRCA1_1.tsv", "T_BRCA1_2.tsv", "T_BRCA1_3.tsv", "T_BRCA1_4.tsv", "T_BRCA1_5.tsv", 
             "T_BRCA2_1.tsv", "T_BRCA2_2.tsv", "T_BRCA2_3.tsv", "T_BRCA2_4.tsv", "T_BRCA2_5.tsv", 
             "T_WT_1.tsv", "T_WT_2.tsv", "T_WT_3.tsv", "T_WT_4.tsv", "T_WT_5.tsv", "T_WT_6.tsv", "T_WT_7.tsv", "T_WT_8.tsv", "T_WT_9.tsv", "T_WT_10.tsv")
sampleName <- c("T_BRCA1_1", "T_BRCA1_2", "T_BRCA1_3", "T_BRCA1_4", "T_BRCA1_5", 
                "T_BRCA2_1", "T_BRCA2_2", "T_BRCA2_3", "T_BRCA2_4", "T_BRCA2_5", 
                "T_WT_1", "T_WT_2", "T_WT_3", "T_WT_4", "T_WT_5", "T_WT_6", "T_WT_7", "T_WT_8", "T_WT_9", "T_WT_10")
Condition <- c('BRCA1', 'BRCA1', 'BRCA1', 'BRCA1', 'BRCA1',
               'BRCA2', 'BRCA2', 'BRCA2', 'BRCA2', 'BRCA2',
               'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT', 'WT')
sampleTable <- data.frame(sample=sampleName, file=fileName, condition=Condition)
files <- file.path(directory, sampleTable$file)
names(files) <- paste0(sampleName)

library(tximport) # tximport_1.14.2

txdb <- read.table("txdb_geneSymbol_onlyProtein.txt", header=TRUE)
txi <- tximport(files, type="kallisto", tx2gene=txdb)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ condition)
dds <- ddsTxi[ rowSums(counts(ddsTxi)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
resLFC <- lfcShrink(dds, type="ashr", res=res)

vsd <- vst(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup = c("condition"), ntop=5000, returnData=TRUE) # Then, Calculate pearson correlation with gene expression with this PC1 value 
library(ggplot2)
library(ggrepel)
ggplot(a, aes(x=PC1, y=PC2)) + geom_text_repel(label=rownames(pca), size=3)

# HRR and mesenchymal
library('DESeq2')
directory <- # Directory in wich Kallisto tsv files exist
fileName<- c("T_BRCA1_1.tsv", "T_BRCA1_2.tsv", "T_BRCA1_3.tsv", "T_BRCA1_4.tsv", "T_BRCA1_5.tsv", 
             "T_BRCA2_1.tsv", "T_BRCA2_2.tsv", "T_BRCA2_3.tsv", "T_BRCA2_4.tsv", "T_BRCA2_5.tsv", 
             "T_WT_1.tsv", "T_WT_2.tsv", "T_WT_3.tsv", "T_WT_4.tsv", "T_WT_5.tsv", "T_WT_6.tsv", "T_WT_7.tsv", "T_WT_8.tsv", "T_WT_9.tsv", "T_WT_10.tsv")
sampleName <- c("T_BRCA1_1", "T_BRCA1_2", "T_BRCA1_3", "T_BRCA1_4", "T_BRCA1_5", 
                "T_BRCA2_1", "T_BRCA2_2", "T_BRCA2_3", "T_BRCA2_4", "T_BRCA2_5", 
                "T_WT_1", "T_WT_2", "T_WT_3", "T_WT_4", "T_WT_5", "T_WT_6", "T_WT_7", "T_WT_8", "T_WT_9", "T_WT_10")
Condition <- c('HR', 'HR', 'HR', 'HR', 'HR',
               'HR', 'HR', 'Mesenchymal', 'HR', 'HR',
               'Mesenchymal', 'HR', 'HR', 'Mesenchymal', 'HR', 'Mesenchymal', 'Mesenchymal', 'HR', 'HR', 'HR')
sampleTable <- data.frame(sample=sampleName, file=fileName, condition=Condition)
files <- file.path(directory, sampleTable$file)
names(files) <- paste0(sampleName)

library(tximport)

txdb <- read.table("txdb_geneSymbol_onlyProtein.txt", header=TRUE)
txi <- tximport(files, type="kallisto", tx2gene=txdb)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ condition)
dds <- ddsTxi[ rowSums(counts(ddsTxi)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, type="ashr", res=res)

vsd <- vst(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup = c("condition"), ntop=5000, returnData=TRUE)
library(ggplot2)
library(ggrepel)
ggplot(a, aes(x=PC1, y=PC2)) + geom_text_repel(label=rownames(pca), size=3)
