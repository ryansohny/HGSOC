library(survival)
library(survminer)
library(dplyr)

# SNUH data
snuhov <- read.csv("SNUH-OV_clinical.csv")
fit <- survfit(Surv(time.overall, status) ~ emt, data = snuhov)
ggsurvplot(fit, conf.int = TRUE, pval = TRUE, palette = c('#000066','#990000'))

covariates <- c("age", "stage", "lymphnode.invasion", "lvsi", "treatment")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time.overall, status)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = snuhov)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2); #coeficient beta
                         HR <-signif(x$coef[2], digits=2); #exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)

tcgaov <- read.csv("TCGA_clinical.csv")
fit <- survfit(Surv(time, status) ~ emt, data = tcgaov)
ggsurvplot(fit, conf.int = TRUE, pval = TRUE, palette = c('#990000', '#000066'))
