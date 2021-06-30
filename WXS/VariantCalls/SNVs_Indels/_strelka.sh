#!/bin/bash

# Input file ==> Samplelist_Strelka2_Manta.txt

#Usage : ~.sh [start] [end] [sample_list]
python=/usr/bin/python

for((i=$2;i<=$3;i++))
do
    sed -n ${i}p $1 > tmp${i}
    normal=$(awk '{print $2}' tmp${i})
    tumor=$(awk '{print $1}' tmp${i})
    rm tmp${i}
mkdir -p ${tumor}_strelka
cd ${tumor}_strelka
$python ~/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam=${normal}.dp.recal.bam \
--tumorBam=${tumor}.dp.recal.bam \
--referenceFasta human_g1k_v37.fasta \
--indelCandidates ${tumor}_manta/MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz \
--exome \
--callRegions Exome-Agilent_V6_Padded_wo-chr_woY_strelka_manta.bed.gz

cd StrelkaSomaticWorkflow
./runWorkflow.py -m local -j 35
cd ../..
done
