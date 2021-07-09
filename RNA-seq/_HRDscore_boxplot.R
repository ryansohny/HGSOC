library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
dat <- read.csv("EMT_HRD.csv", row.names=1)

png("HRDscore_byCluster.png")
p <- ggboxplot(dat, x="Cluster", y="HRD", color="Cluster", xlab="Cluster", ylab="HRD")
print(p)
dev.off()

png("HRDscore_byCluster.png")


pdf("scarHRD_violoinplot.pdf")
p1 <- ggplot(dat, aes(x=Cluster, y=HRD, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p2 <- ggplot(dat, aes(x=Cluster, y=Telomeric.AI, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p3 <- ggplot(dat, aes(x=Cluster, y=LST, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()
p4 <- ggplot(dat, aes(x=Cluster, y=HRD.sum, fill=Cluster)) + geom_violin(trim=TRUE) + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7) + scale_fill_manual(values=c('#000066','#990000'))  + theme_classic()

ggarrange(p1, p2, p3, p4, labels=c("A", "B", "C", "D"), ncol=2, nrow=2)
dev.off()

#dat %>% wilcox_test(Stromal ~ Cluster, exact=TRUE, detailed=TRUE, p.adjust.method="bonferroni")
