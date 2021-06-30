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

/ruby/Tools/manta/manta-1.5.0/bin/configManta.py \
--normalBam /clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/${normal}/${normal}.dp.recal.bam \
--tumorBam /clinix1/Analysis/mongol/phenomata/09.HGSOC/01.Alignment/${tumor}/${tumor}.dp.recal.bam \
--referenceFasta /clinix1/Users/uugi0620/project/chul/wgs/v37/human_g1k_v37.fasta \
--exome \
--callRegions /clinix1/Analysis/mongol/phenomata/09.HGSOC/Database/Exome-Agilent_V6_Padded_wo-chr_woY_strelka_manta.bed.gz

cd MantaWorkflow
./runWorkflow.py -m local -j 35
cd ../..

done
