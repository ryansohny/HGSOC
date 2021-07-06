# This code needs to be more sophisticated

library(data.table)
dat <- read.csv("/clinix1/Analysis/mongol/phenomata/09.HGSOC/Genes_Revision/SomaticMutation_Mutual_exclusiveness.csv", row.names=1)
dat2 <- dat[2:101,] # Except for TP53 somatic mutation
sample_list <- list("T_BRCA1_1", "T_BRCA1_2", "T_BRCA1_3", "T_BRCA1_4", "T_BRCA1_5", "T_BRCA2_1", "T_BRCA2_2", "T_BRCA2_3", "T_BRCA2_4", "T_BRCA2_5", "T_WT_1", "T_WT_2", "T_WT_3", "T_WT_4", "T_WT_5", "T_WT_8", "T_WT_9", "T_WT_10") # T_WT_6 and T_WT_7 ==> NO somatic SNVs/Indels
comparison <- list()

# Total of 153 pairwise Fisher's Exact test
for(i in seq(1,17)){
        for(j in seq(i+1,18)){
                p_fisher <- fisher.test(table(as.factor(dat2[,sample_list[[i]]]), as.factor(dat2[,sample_list[[j]]])))
                comparison <- append(comparison, paste(paste(sample_list[i], sample_list[j], sep=":"), as.character(p_fisher$p.value), sep="\t"))
        }
}

# Check for Mutual Exclusiveness
for(i in comparison){
        cat(paste(i, "\n", sep=""))
}
