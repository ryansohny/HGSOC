library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
dat <- read.csv("HRD_score.csv", row.names=1)

dat %>% wilcox_test(HRD ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic       p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>    <dbl>     <dbl>
#1     14.0 HRD   A      B         15     5      72.5 0.00256     7.00      18.0
dat %>% wilcox_test(Telomeric.AI ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic      p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl>    <dbl>     <dbl>
#1     17.0 Telo… A      B         15     5        61 0.0444     2.00      24.0
dat %>% wilcox_test(LST ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic     p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>     <dbl>
#1     21.0 LST   A      B         15     5        66 0.014     2.00      27.0
dat %>% wilcox_test(HRD.sum ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
## A tibble: 1 x 12
#  estimate .y.   group1 group2    n1    n2 statistic      p conf.low conf.high
#*    <dbl> <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl>    <dbl>     <dbl>
#1     51.0 HRD.… A      B         15     5        67 0.0113     11.0      70.0

pdf("scarHRD_violoinplot.pdf")
p1 <- ggplot(dat, aes(x=Cluster, y=HRD, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p2 <- ggplot(dat, aes(x=Cluster, y=Telomeric.AI, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p3 <- ggplot(dat, aes(x=Cluster, y=LST, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p4 <- ggplot(dat, aes(x=Cluster, y=HRD.sum, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()

ggarrange(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
dev.off()
