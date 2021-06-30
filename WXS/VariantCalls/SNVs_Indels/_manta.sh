#!/bin/bash

# Input file ==> Samplelist_Strelka2_Manta.txt
for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        normal=$(awk '{print $2}' tmp${i})
        tumor=$(awk '{print $1}' tmp${i})
        rm tmp${i}

mkdir -p ${tumor}_manta
cd ${tumor}_manta

~/manta-1.5.0/bin/configManta.py \
--normalBam ${normal}.dp.recal.bam \
--tumorBam ${tumor}.dp.recal.bam \
--referenceFasta human_g1k_v37.fasta \
--exome \
--callRegions Exome-Agilent_V6_Padded_wo-chr_woY_strelka_manta.bed.gz

cd MantaWorkflow
./runWorkflow.py -m local -j 35
cd ../..

done
