source activate sequenza
reference=human_g1k_v37.fasta
gcwiggle=v37.gc50Base.wig.gz
# Input File ==> samplelist_sequenza.txt

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        normal=$(awk '{print $3}' tmp${i})
        tumor=$(awk '{print $2}' tmp${i})
        rm -rf tmp${i}

mkdir -p $sample
echo 'Generation of seqz file : ' ${sample}
sequenza-utils bam2seqz \
-n ${normal} \
-t ${tumor} \
--fasta $reference \
-gc $gcwiggle \
-o ${sample}/${sample}.seqz.gz

echo 'Post-process by binning the original seqz file : ' ${sample}
sequenza-utils seqz_binning \
--seqz ${sample}/${sample}.seqz.gz \
-w 50 \
-o ${sample}/${sample}.w50.seqz.gz
done
