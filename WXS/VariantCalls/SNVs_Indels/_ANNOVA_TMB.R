# TMB ANOVA test
df <- data.frame(group=c("BRCA1", "BRCA1", "BRCA1", "BRCA1", "BRCA1", "BRCA2", "BRCA2", "BRCA2", "BRCA2", "BRCA2", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT"), TMB=c(2.75247524752475,	5.89108910891089,	5.47524752475248,	2.95049504950495,	6.85148514851485,	6.07920792079208,	4.43564356435644,	0.821782178217822,	3.7029702970297,	3.46534653465347,	2.84158415841584,	5.38613861386139,	1.88118811881188,	2.25742574257426,	2.81188118811881,	0.821782178217822,	0.693069306930693,	3.5049504950495,	4.72277227722772,	6.53465346534653))
levels(df$group)
library("ggpubr")
ggboxplot(df, x = "group", y = "TMB", 
          color = "group", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("BRCA1", "BRCA2", "WT"),
          ylab = "TMB", xlab = "Group")
ggline(df, x = "group", y = "TMB", 
       add = c("mean_se", "jitter"), 
       order = c("BRCA1", "BRCA2", "WT"),
       ylab = "TMB", xlab = "Group")

res.aov <- aov(TMB ~ group, data = df)
summary(res.aov)
