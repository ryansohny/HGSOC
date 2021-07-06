#!/bin/bash

source activate gatk4

reference=human_g1k_v37.fasta
region=Exome-Agilent_V6_Padded_wo-chr_woY.bed
Mills=Mills_and_1000G_gold_standard.indels.b37.vcf
KG=1000G_phase1.indels.b37.vcf
dbsnp=dbsnp_138.b37.vcf
threads=

# Step 1 : Run Mutect2 in tumor-only mode for each normal sample
for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        rm tmp${i}
mkdir -p ${sample}
$gatk4 --java-options "-Xmx100G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
HaplotypeCaller \
-R $reference \
-I ${sample}.dp.recal.bam \
-ERC GVCF \
-L $region \
-O ./${sample}/${sample}.g.vcf \
--native-pair-hmm-threads 10 \
-pairHMM FASTEST_AVAILABLE
done

# Combine multiple GVCFs from HaplotypeCaller
$gatk4 --java-options "-Xmx100G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
CombineGVCFs \
--dbsnp $dbsnp \
-L $region \
-R $reference \
-V N_BRCA1_1.g.vcf \
-V N_BRCA1_2.g.vcf \
-V N_BRCA1_3.g.vcf \
-V N_BRCA1_4.g.vcf \
-V N_BRCA1_5.g.vcf \
-V N_BRCA2_1.g.vcf \
-V N_BRCA2_2.g.vcf \
-V N_BRCA2_3.g.vcf \
-V N_BRCA2_4.g.vcf \
-V N_BRCA2_5.g.vcf \
-V N_WT_1.g.vcf \
-V N_WT_2.g.vcf \
-V N_WT_3.g.vcf \
-V N_WT_4.g.vcf \
-V N_WT_5.g.vcf \
-V N_WT_6.g.vcf \
-V N_WT_7.g.vcf \
-V N_WT_8.g.vcf \
-V N_WT_9.g.vcf \
-V N_WT_10.g.vcf \
-O HGOC_HC.variants.raw.vcf

# Genotype GVCFs
$gatk4 --java-options "-Xmx100G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
GenotypeGVCFs \
--dbsnp $dbsnp \
-L $region \
-R $reference \
-V HGOC_HC.variants.raw.g.vcf \
-O HGOC_HC.variants.raw.vcf

# Select variant : SNPs
$gatk4 --java-options "-Xmx200G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
SelectVariants \
--select-type-to-include SNP \
-L $region \
-R $reference \
-V HGOC_HC.variants.raw.vcf \
-O HGOC_HC.variants.raw.snps.vcf

# Select variant : INDELs
$gatk4 --java-options "-Xmx200G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
SelectVariants \
--select-type-to-include INDEL \
-L $region \
-R $reference \
-V HGOC_HC.variants.raw.vcf \
-O HGOC_HC.variants.raw.indels.vcf

# Variant Filtration : SNPs
$gatk4 --java-options "-Xmx200G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
VariantFiltration \
-L $region \
-R $reference \
-V HGOC_HC.variants.raw.snps.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "GGMI_SNP_FILTER" \
-O HGOC_HC.variants.filtered.snps.vcf

# Variant Filtration : INDELs
$gatk4 --java-options "-Xmx200G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
VariantFiltration \
-L $region \
-R $reference \
-V HGOC_HC.variants.raw.indels.vcf \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "GGMI_INDEL_FILTER" \
--cluster-size 2 \
--cluster-window-size 5 \
-O HGOC_HC.variants.filtered.indels.vcf

# Select Variant SNPS and INDELs
$gatk4 --java-options "-Xmx200G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
SelectVariants \
-L $region \
-R $reference \
-V HGOC_HC.variants.filtered.snps.vcf \
--select-type-to-include SNP \
-O HGOC_HC.variants.filtered.snps.only.vcf

$gatk4 --java-options "-Xmx200G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
SelectVariants \
-L $region \
-R $reference \
-V HGOC_HC.variants.filtered.indels.vcf \
--select-type-to-include INDEL \
-O HGOC_HC.variants.filtered.indels.only.vcf

# Combine Filtered SNP & INDEL VCFs
$gatk4 --java-options "-Xmx200G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Djava.io.tmpdir=tmp" \
MergeVcfs \
-R $reference \
-I HGOC_HC.variants.filtered.snps.only.vcf \
-I HGOC_HC.variants.filtered.indels.only.vcf \
-O HGOC_HC.variants.filtered.combined.vcf

# Run Annovar
annovar_folder="annovar/"
annovar_db_folder="annovar/humandb/"

mkdir -p Annotation

perl ${annovar_folder}table_annovar.pl Annotation/HGOC_HC.annovar ${annovar_db_folder} --buildver hg19 --protocol refGene,knownGene,ensGene,clinvar_20200316,1000g2015aug_all,gnomad211_genome,gnomad211_exome,avsnp150,exac03nontcga,icgc21,cosmic70,nard_v1-1,dbnsfp35c,wgEncodeBroadHmmGm12878HMM,tfbsConsSites --operation g,g,g,f,f,f,f,f,f,f,f,f,f,f,f --remove --outfile Annotation/alltranscript_HGOC_HC -otherinfo
