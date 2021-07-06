source activate cnvkit

# Input file ==> samplelist_cnvkit.txt
for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        rm tmp${i}
cnvkit.py segmetrics -s ${sample}.dp.cn{s,r} --ci -o ${sample}.dp.segmetrics.cns
cnvkit.py call ${sample}.dp.segmetrics.cns --filter ci -o ${sample}.somatic.segmetrics.call.cns
done
