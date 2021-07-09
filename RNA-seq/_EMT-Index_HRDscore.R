library(ggpubr)
library(dplyr)
library(rstatix)

dat <- read.table("XCell_results.txt", header=TRUE, row.names=1)

## Mann-Whitney U test of each cell type enrichment result between cluster A and cluster B ##
dat %>% wilcox_test(Stromal ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")



# Draw Scatterplot
png("Immune_Stromal_scatterplot.png")
ggscatter(dat, x="Immune", y="Stromal", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="spearman", xlab="Immune score", ylab="Stromal score")
dev.off()
