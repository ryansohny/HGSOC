mkdir -p SomaticCNA

cnvkit.py batch \
T_*.dp.recal.bam \
--normal N_*dp.recal.bam \
-t Exome-Agilent_V6_wochr_woY.bed \
-f human_g1k_v37_Autosome_X.fasta \ # Reference Fasta without Y chromosome
-g access-5k-mappable.grch37_woY.bed \
--output-reference ./SomaticCNA/HGSOC.cnn \
--output-dir SomaticCNA/ \
--scatter \
--diagram \
-p 39
