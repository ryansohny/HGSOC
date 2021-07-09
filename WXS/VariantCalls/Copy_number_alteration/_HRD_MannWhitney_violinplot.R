library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
dat <- read.csv("HRD_score.csv", row.names=1)

dat %>% wilcox_test(HRD ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
dat %>% wilcox_test(Telomeric.AI ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
dat %>% wilcox_test(LST ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
dat %>% wilcox_test(HRD.sum ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")

pdf("scarHRD_violoinplot.pdf")
p1 <- ggplot(dat, aes(x=Cluster, y=HRD, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p2 <- ggplot(dat, aes(x=Cluster, y=Telomeric.AI, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p3 <- ggplot(dat, aes(x=Cluster, y=LST, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p4 <- ggplot(dat, aes(x=Cluster, y=HRD.sum, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()

ggarrange(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
dev.off()
