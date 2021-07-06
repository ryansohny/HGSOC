source activate gistic2

basedir="GISTIC2/"
mkdir -p $basedir

echo --- running GISTIC ---
segfile=HGSOC_CNVkit.seg

refgenefile=hg19.mat
gistic2="gistic2 bin executable"

$gistic2 \
-b $basedir \
-seg $segfile \
-refgene $refgenefile \
-ta 0.1 \
-td 0.1 \
-js 100 \
-genegistic 1 \
-conf 0.90 \
-savegene 1
