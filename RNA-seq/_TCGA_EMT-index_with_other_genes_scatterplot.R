library(ggpubr)
library(dplyr)
library(rstatix)
library(data.table)

dat <- read.table("/clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/RNA-seq/Analysis_new/GDC_ovarian_cancer/Z_Analysis/TPM/Z_GeneID/Fig5_Pearson_EMT-index_Genes_TCGA-OV_TPM_GeneLevel_geneID_log2.txt", header=TRUE, row.names=1)
t_dat <- transpose(dat)

colnames(t_dat) <- rownames(dat)
rownames(t_dat) <- colnames(dat)

### Draw Scatterplot
# Figure 5A
pdf("Scatterplot_EMT-index_Figure5A.pdf")
p1 <- ggscatter(t_dat, x="EMT_index", y="CDH1", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="ZEB1")
p2 <- ggscatter(t_dat, x="EMT_index", y="CDH2", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="CDH2")
p3 <- ggscatter(t_dat, x="EMT_index", y="VIM", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="VIM")
p4 <- ggscatter(t_dat, x="EMT_index", y="TGFB1", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="TGFB1")
ggarrange(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
dev.off()


pdf("Scatterplot_EMT-index_EMT5TFs.pdf", 11, 7)
p1 <- ggscatter(t_dat, x="EMT_index", y="TWIST1", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="TWIST1")
p2 <- ggscatter(t_dat, x="EMT_index", y="SNAI1", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="SNAI1")
p3 <- ggscatter(t_dat, x="EMT_index", y="SNAI2", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="SNAI2")
p4 <- ggscatter(t_dat, x="EMT_index", y="ZEB1", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="ZEB1")
p5 <- ggscatter(t_dat, x="EMT_index", y="ZEB2", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="pearson", xlab="EMT index", ylab="ZEB2")
ggarrange(p1, p2, p3, p4, p5, labels=c("A", "B", "C", "D", "E"), ncol=3, nrow=2)
dev.off()

