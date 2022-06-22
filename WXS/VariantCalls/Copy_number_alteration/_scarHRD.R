library(scarHRD)

sample_list <- list("T_BRCA1_1", "T_BRCA1_2", "T_BRCA1_3", "T_BRCA1_4", "T_BRCA1_5", "T_BRCA2_1", "T_BRCA2_2", "T_BRCA2_3", "T_BRCA2_4", "T_BRCA2_5", "T_WT_1", "T_WT_2", "T_WT_3", "T_WT_4", "T_WT_5", "T_WT_6", "T_WT_7", "T_WT_8", "T_WT_9", "T_WT_10")

dir <- "/clinix1/Analysis/mongol/phenomata/09.HGSOC/02.VariantCalls/Sequenza"

for(sample in sample_list) {
        datadir <- paste(dir, sample, sample, sep='/')
        data.file <- paste(datadir, ".w50.seqz.gz", sep='')
        scar_score(data.file, reference="grch37", chr.in.names=FALSE, seqz=TRUE)
        }
