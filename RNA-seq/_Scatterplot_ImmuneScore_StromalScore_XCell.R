library(ggpubr)
dat <- read.table("ImmuneScore_StromalScore_XCell.txt")

# Wilcox.test


# Draw Scatterplot
png("Immune_Stromal_scatterplot.png")
ggscatter(dat, x="Immune", y="Stromal", add="reg.line", conf.int=TRUE, cor.coef=TRUE, cor.method="spearman", xlab="Immune score", ylab="Stromal score")
dev.off()
