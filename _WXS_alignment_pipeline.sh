#!/bin/bash

bwa=bwa exec # version 0.7.17
samtools=samtools exec # version 1.9
gatk4=gatk exec # version 4.1.4.1
reference=human_g1k_v37.fasta
region=Exome-Agilent_V6_Padded_wo-chr_woY.bed
Mills=Mills_and_1000G_gold_standard.indels.b37.vcf
indels=1000G_phase1.indels.b37.vcf
dbsnp=dbsnp_138.b37.vcf
threads=

# input file ==> Sample_WXS.txt
for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        read1=$(awk '{print $2}' tmp${i})
        read2=$(awk '{print $3}' tmp${i})
        rm tmp${i}
mkdir -p ${sample}
line=`zcat $fq1 | head -1`
tagID=`echo "$line" | awk -F":" '{print $1,":",$2,":",$3}' | sed 's/@//g' | sed 's/ //g'`
platform="NextSeq"
library=$sample
rgline="@RG\tPL:$platform\tID:$tagID\tSM:$sample\tLB:$library"

# BWA Alignment and Sorting 
$bwa mem -M -t $threads -R "$rgline" $reference $read1 $read2 \
| $samtools sort -@ $threads -O BAM -o ${sample}/${sample}.bam -

# Mark Duplicates
mkdir -p Metrics
$gatk4 MarkDuplicates \
-I ./${sample}/${sample}.bam \
-M ./Metrics/${sample}.metrics \
-ASSUME_SORT_ORDER coordinate \
-REMOVE_DUPLICATES true \
-CREATE_INDEX true \
-O ./${sample}/${sample}.dp.bam

# Base Recalibrator
$gatk4 BaseRecalibrator \
-I ./${sample}/${sample}.dp.bam \
-R ${reference} \
--known-sites ${dbsnp} \
--known-sites ${indels} \
--known-sites ${Mills} \
-L ${region} \
-O ./${sample}/${sample}.tabl

# Apply Base quality score recalibration
$gatk4 ApplyBQSR \
-bqsr ./${sample}/${sample}.tabl \
-I ./${sample}/${sample}.dp.bam \
-O ./${sample}/${sample}.dp.recal.bam \
-L ${region}

# Depth of Coverage
mkdir -p Depth
$samtools depth \
${sample}/${sample}.dp.recal.bam \
-q 20 \
-Q 20 \
-b ${region} \
|  awk '{sum+=$3} END { print sum/NR}' \
> Depth/${sample}.depth

rm ${sample}/${sample}.dp.bam
rm ${sample}/${sample}.bam

done
