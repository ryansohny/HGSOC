#python3
from glob import glob
import subprocess as sp
a = sorted(glob('*.targetCoverage.out.sample_summary'))
rfh = open("WES_QC_metrics.txt", 'w')
rfh.write("Sample ID\tTotal number of covered bases\tMedian coverage per base\t3rd quantile coverage per base\t1st quantile coverage per base\tPercentage of targeted bases with coverage>=10\tPercentage of targeted bases with coverage>=20\tPercentage of targeted bases with coverage>=30\n")
rfh.flush()
for i in a:
        id = i.split('.')[0]
        covered = str(int(sp.getoutput("awk '$2 > 0 { count++} END { print count }' " + id + ".targetCoverage.out"))-1)
        dfh = open(i, 'r')
        dfh.readline()
        line = dfh.readline().strip().split('\t')
        rfh.write(id + '\t' + covered + '\t' + line[4] + '\t' + line[3] + '\t' + line[5] + '\t' + line[7] + '\t' + line[8] + '\t' + line[9] + '\n')
        rfh.flush()
        dfh.close()
