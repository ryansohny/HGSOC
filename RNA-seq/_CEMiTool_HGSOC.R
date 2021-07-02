library(CEMiTool)

input <- # VST applied HGSOC RNA-seq data table
expr <- read.table(input, row.names = 1, sep="\t", header=T)
sample_annot <- read.csv("CEMiTool_annotation.csv")
int_string <- read.csv("ARCHS4_Coexpression_interaction.csv")
hallmark_gmt <- read_gmt("h.all.v7.4.symbols.gmt") # http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/h.all.v7.4.symbols.gmt


# Module detection
cem <- cemitool(expr=expr, annot=sample_annot, gmt=hallmark_gmt, cor_method="spearman", interactions=int_string, merge_similar=TRUE, network_type="signed", tom_type="signed", rank_method="mean",class_column="cluster", gsea_max_siz=2000, verbose=TRUE)

write_files(cem, directory="./Tables", force=TRUE)

# Plotting

cem <- mod_gsea(cem)
show_plot(cem, "gsea")
dev.off()

cem <- mod_ora(cem, go_gmt)
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
pdf("ORA.pdf")
plots
dev.off()

library(ggplot2)
interactions_data(cem) <- int_string
cem <- plot_interactions(cem)
plots <- show_plot(cem, "interaction")
pdf("interaction_ARCHS4.pdf")
plots
dev.off()
