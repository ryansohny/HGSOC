#!/bin/bash

kallisto=<kallisto exec> # version 0.46.1
index=GRCh37_rna.idx

# input_file ==> Sample_RNA.txt
for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        read1=$(awk '{print $2}' tmp${i})
        read2=$(awk '{print $3}' tmp${i})
        rm tmp${i}
$kallisto quant --index=${index} --threads 40 --output-dir=${sample} ${read1} ${read2}
done
