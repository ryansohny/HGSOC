import gzip as gz
import sys
dbf = gz.open("GRCh37_latest_genomic.gtf.gz", 'rb') # http://mirror.ufs.ac.za/datasets/ncbi/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gtf.gz

dfh = open(sys.argv[1], 'r')
rfh = open(sys.argv[1][:-4] + '_geneSymbol.txt', 'w')

dbf.readline()
dbf.readline()
dbf.readline()
dbf.readline()
dbf.readline()
db = {}
for i in dbf:
	line = i.decode().strip().split('\t')
	try:
		if line[8].split(';')[1].split('"')[0].strip() == 'transcript_id':
			if line[8].split(';')[1].split('"')[1] not in db.keys():
				db[line[8].split(';')[1].split('"')[1]] = line[8].split(';')[0].split('"')[1]
	except IndexError:
		pass
line = dfh.readline()
rfh.write(line)
for i in dfh:
	line = i.strip().split('\t')
	rfh.write(db[line[0]] + '\t' + '\t'.join(line[1:]) + '\n')
	rfh.flush()
