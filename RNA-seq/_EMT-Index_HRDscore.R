library(ggpubr)
library(dplyr)
library(rstatix)

dat <- read.csv("EMT_HRD.csv", header=TRUE, row.names=1)

# Draw Scatterplot
pdf("EMT-index_HRD_scatterplot.pdf")
ggscatter(dat, x="EMT.index", y="HRD", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="spearman", xlab="EMT index", ylab="HRD score")
dev.off()

