import sys
dfh = open(sys.argv[1], 'r')
rfh = open(sys.argv[1].rstrip('geneSymbol.txt') + 'GeneLevel_geneSymbol.txt', 'w')

rfh.write(dfh.readline())
exp_dic = {}
for i in dfh:
	line = i.strip().split('\t')
	if line[0] not in exp_dic:
		exp_dic[line[0]] = '\t'.join(line[1:])
	else:
		new = []
		for j in range(0, len(line[1:])):
			new.append(str(float(exp_dic[line[0]].split('\t')[j]) + float(line[j+1])))
		exp_dic[line[0]] = '\t'.join(new)

for key in exp_dic:
	rfh.write(key + '\t' + exp_dic[key] + '\n')
	rfh.flush()

dfh.close()
rfh.close()
