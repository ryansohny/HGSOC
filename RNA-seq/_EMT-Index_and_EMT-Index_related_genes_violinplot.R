library(ggpubr)
library(dplyr)
library(rstatix)
library(data.table)
library(ggarrange)

dat <- read.table("EMT_Index_HGSOC.txt", header=TRUE, row.names=1)
dat %>% wilcox_test(EMT_Index ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic       p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>    <dbl>     <dbl>
#1    -7.47 EMT_… A      B         15     5         0 1.29e-4    -9.70     -4.50

pdf("Violinplot_EMT-Index_Figure3F.pdf")
ggplot(dat, aes(x=Cluster, y=EMT_Index, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
dev.off()

dat <- read.table("CDH1_VIM_TGFB1_HGOC_TPM_new_GeneLevel_geneSymbol_onlyTumor_log2.txt", header=TRUE, row.names=1)
t_dat <- transpose(dat)

colnames(t_dat) <- rownames(dat)
rownames(t_dat) <- colnames(dat)

t_dat[,"TGFB1"] <- as.numeric(t_dat[,"TGFB1"])
t_dat[,"CDH1"] <- as.numeric(t_dat[,"CDH1"])
t_dat[,"VIM"] <- as.numeric(t_dat[,"VIM"])
t_dat[,"Cluster"] <- as.factor(t_dat[,"Cluster"])

t_dat %>% wilcox_test(CDH1 ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic       p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>    <dbl>     <dbl>
#1     2.86 CDH1  A      B         15     5        73 5.16e-4     1.41      4.08
## … with 2 more variables: method <chr>, alternative <chr>

t_dat %>% wilcox_test(VIM ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic       p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>    <dbl>     <dbl>
#1    -3.01 VIM   A      B         15     5         0 1.29e-4    -3.97     -2.01
## … with 2 more variables: method <chr>, alternative <chr>

t_dat %>% wilcox_test(TGFB1 ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic       p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>    <dbl>     <dbl>
#1    -1.14 TGFB1 A      B         15     5         3 9.03e-4    -1.78    -0.631
## … with 2 more variables: method <chr>, alternative <chr>

pdf("Violinplot_CDH1-VIM-TGFB1_Figure3F.pdf")
p1 <- ggplot(t_dat, aes(x=Cluster, y=CDH1, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p2 <- ggplot(t_dat, aes(x=Cluster, y=VIM, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p3 <- ggplot(t_dat, aes(x=Cluster, y=TGFB1, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()

ggarrange(p1, p2, p3, labels=c("A", "B", "C"), ncol=2, nrow=2)
dev.off()
