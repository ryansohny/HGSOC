#!/bin/bash

gatk4=gatk exec # version 4.1.4.1

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        read1=$(awk '{print $2}' tmp${i})
        read2=$(awk '{print $3}' tmp${i})
        rm tmp${i}
$gatk4 \
-T DepthOfCoverage \
-R human_g1k_v37.fasta \
-I ${sample}.dp.recal.bam \
-L Exome-Agilent_V6_Padded_wo-chr_woY.bed \
-o ${sample}.targetCoverage.out \
-baseCounts \
-dt None \
--minMappingQuality 20 \
--minBaseQuality 20 \
--logging_level ERROR \
--summaryCoverageThreshold 1 \
--summaryCoverageThreshold 10 \
--summaryCoverageThreshold 20 \
--summaryCoverageThreshold 30
done
