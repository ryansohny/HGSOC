library(ggpubr)
library(dplyr)
library(rstatix)
library(data.table)
library(ggarrange)

dat <- read.table("CTNNB1_TCF4_ZEB1_HGOC_TPM_new_GeneLevel_geneSymbol_onlyTumor_log2.txt", header=TRUE, row.names="ID")
t_dat <- transpose(dat)
colnames(t_dat) <- rownames(dat)
rownames(t_dat) <- colnames(dat)

t_dat[,"CTNNB1"] <- as.numeric(t_dat[,"CTNNB1"])
t_dat[,"TCF4"] <- as.numeric(t_dat[,"TCF4"])
t_dat[,"ZEB1"] <- as.numeric(t_dat[,"ZEB1"])

dat %>% wilcox_test(CTNNB1 ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic      p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl>    <dbl>     <dbl>
#1   -0.532 CTNN… A      B         15     5        11 0.0193   -0.862   -0.0869
dat %>% wilcox_test(TCF4 ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic       p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>    <dbl>     <dbl>
#1    -1.69 TCF4  A      B         15     5         0 1.29e-4    -2.14     -1.19
dat %>% wilcox_test(ZEB1 ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic       p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>    <dbl>     <dbl>
#1    -2.11 ZEB1  A      B         15     5         0 1.29e-4    -2.87     -1.35

pdf("Violinplot_CTNNB1-TCF4-ZEB1_FigureS11.pdf",16,7)

p1 <- ggplot(t_dat, aes(x=Cluster, y=CTNNB1, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p2 <- ggplot(t_dat, aes(x=Cluster, y=TCF4, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p3 <- ggplot(t_dat, aes(x=Cluster, y=ZEB1, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()

ggarrange(p1, p2, p3, labels=c("A", "B", "C"), ncol=3, nrow=1)
dev.off()
