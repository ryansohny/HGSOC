#!python2

import sys, operator
import gzip as gz
from subprocess import call
from collections import OrderedDict
dfh = gz.open(sys.argv[1], 'rb') # MAF file of all tumors are merged into one gzipped file and used as an input
rfh1 = open(sys.argv[1][:-7] + '_matrix.txt', 'w') # 0: no mut, 1: yes mut
rfh2 = open(sys.argv[1][:-7] + '_matrix_detail.txt', 'w') # detailed mutational classfication, TMB (filter criteria needed), mutational spectrum...
rfh3 = open(sys.argv[1][:-7] + '_TMB.txt', 'w') # Tumor Mutational burden
rfh4 = open(sys.argv[1][:-7] + '_sig.txt', 'w') # Signature
rfh5 = open(sys.argv[1][:-7] + '_titv.txt', 'w') # TiTv
dfh.readline()
# full_db = ["3'UTR", "5'Flank", "5'UTR", 'De_novo_Start_InFrame', 'De_novo_Start_OutOfFrame', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'IGR', 'In_Frame_Del', 'In_Frame_Ins', 'Intron', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'RNA', 'Silent', 'Splice_Site', 'Start_Codon_Del', 'Start_Codon_Ins', 'Start_Codon_SNP', 'Stop_Codon_Del', 'Stop_Codon_Ins', 'lincRNA']
exonic_db = {'Missense_Mutation':'1', 'Nonsense_Mutation':'2', 'Frame_Shift_Del':'3', 'Frame_Shift_Ins':'4', 'In_Frame_Del':'5', 'In_Frame_Ins':'6', 'Splice_Site':'7'}
from collections import OrderedDict
# line[0] ==> Gene Name
# line[8] ==> Variant_Classification
# line[15] ==> Tumor_sample_Barcode
t_mutc = dict() # Total mutation counts (not TMB, this is for exonic)
check_dup = list()
line = dfh.readline().strip('\n').split('\t')
sample = line[15]
if line[8] in exonic_db.iterkeys():
	t_mutc[line[0]] = 1
	check_dup.append(line[0])
line = dfh.readline().strip('\n').split('\t')
while line != ['']:
        if line[15] == sample:
		if line[8] in exonic_db.iterkeys():
			if line[0] not in check_dup:
				check_dup.append(line[0])
				try:
					t_mutc[line[0]] += 1
				except KeyError:
					t_mutc[line[0]] = 1
	else:
		check_dup = list()
		if line[8] in exonic_db.iterkeys():
			check_dup.append(line[0])
			try:
				t_mutc[line[0]] += 1
			except KeyError:
				t_mutc[line[0]] = 1
		sample = line[15]
	line = dfh.readline().strip('\n').split('\t')

sorted_t_mutc = sorted(t_mutc.items(), key=operator.itemgetter(1), reverse=True)

dict_stm = OrderedDict(sorted_t_mutc) # dictionary of sorted_t_mutc
dict_sm = OrderedDict(map(lambda x: (x, '0'), map(lambda y: y[0], sorted_t_mutc))) # dictionary of mutation presense of each sample
dict_dsm = OrderedDict(map(lambda x: (x, '0'), map(lambda y: y[0], sorted_t_mutc))) # dictionary of "detailed" mutation presence of each sample
dfh.seek(0)
dfh.readline()
def od_to_tuple(ordered_dict):
	return tuple(map(lambda x: str(ordered_dict[x]), ordered_dict.iterkeys()))
unzipped_mutinfo = list() # later zip() into integrated list to form comutation matrix
unzipped_dmutinfo = list() # detailed mut info
tmb = dict()
#### First line ####
sigdb = OrderedDict([("ACA_CA",0), ("ACC_CA",0), ("ACG_CA",0), ("ACT_CA",0), ("CCA_CA",0), ("CCC_CA",0), ("CCG_CA",0), ("CCT_CA",0), ("GCA_CA",0), ("GCC_CA",0), ("GCG_CA",0), ("GCT_CA",0), ("TCA_CA",0), ("TCC_CA",0), ("TCG_CA",0), ("TCT_CA",0), ("ACA_CG",0), ("ACC_CG",0), ("ACG_CG",0), ("ACT_CG",0), ("CCA_CG",0), ("CCC_CG",0), ("CCG_CG",0), ("CCT_CG",0), ("GCA_CG",0), ("GCC_CG",0), ("GCG_CG",0), ("GCT_CG",0), ("TCA_CG",0), ("TCC_CG",0), ("TCG_CG",0), ("TCT_CG",0), ("ACA_CT",0), ("ACC_CT",0), ("ACG_CT",0), ("ACT_CT",0), ("CCA_CT",0), ("CCC_CT",0), ("CCG_CT",0), ("CCT_CT",0), ("GCA_CT",0), ("GCC_CT",0), ("GCG_CT",0), ("GCT_CT",0), ("TCA_CT",0), ("TCC_CT",0), ("TCG_CT",0), ("TCT_CT",0), ("ATA_TA",0), ("ATC_TA",0), ("ATG_TA",0), ("ATT_TA",0), ("CTA_TA",0), ("CTC_TA",0), ("CTG_TA",0), ("CTT_TA",0), ("GTA_TA",0), ("GTC_TA",0), ("GTG_TA",0), ("GTT_TA",0), ("TTA_TA",0), ("TTC_TA",0), ("TTG_TA",0), ("TTT_TA",0), ("ATA_TC",0), ("ATC_TC",0), ("ATG_TC",0), ("ATT_TC",0), ("CTA_TC",0), ("CTC_TC",0), ("CTG_TC",0), ("CTT_TC",0), ("GTA_TC",0), ("GTC_TC",0), ("GTG_TC",0), ("GTT_TC",0), ("TTA_TC",0), ("TTC_TC",0), ("TTG_TC",0), ("TTT_TC",0), ("ATA_TG",0), ("ATC_TG",0), ("ATG_TG",0), ("ATT_TG",0), ("CTA_TG",0), ("CTC_TG",0), ("CTG_TG",0), ("CTT_TG",0), ("GTA_TG",0), ("GTC_TG",0), ("GTG_TG",0), ("GTT_TG",0), ("TTA_TG",0), ("TTC_TG",0), ("TTG_TG",0), ("TTT_TG", 0)])
titvdb = OrderedDict([("CT",0),("TC",0),("TG",0), ("TA",0), ("CG",0), ("CA",0)])
#line[65] ==> ref context, Format ==> 10 nucleotides [ref allele] 10 nucleotides
# rfh4 ==> sig
# rfh5 ==> titv
# return complementary 5'-XYZ-3' of input trinucleotide if a Y nucleotide is either 'A' or 'G' (purine)
rfh4.write('Sample\t' + '\t'.join(sigdb.iterkeys()) + '\n')
rfh5.write('Sample\tG>A:C>T(Ti)\tA>G:T>C(Ti)\tA>C:T>G(Tv)\tA>T:T>A(Tv)\tG>C:C>G(Tv)\tG>T:C>A(Tv)\n')
def comp(string):
        if string[1] == 'A' or string[1] =='G':
                return string.replace("A",'a').replace("T",'t').replace("G",'g').replace('C','c').replace('a','T').replace('t','A').replace('g','C').replace('c','G')[::-1]
        else:
                return string
#rfh.write('Category\tMutations/Mb\t' + '\t'.join(totalsig.keys()) + '\n')
#rfh2.write('Category\tMutations/Mb\tG.C>A.T(Ti)\tA.T>G.C(Ti)\tA.T>C.G(Tv)\tA.T>T.A(Tv)\tG.C>C.G(Tv)\tG.C>T.A(Tv)\n')
sigsub = {"CA":"CA", "CG":"CG", "CT":"CT", "TA":"TA", "TC":"TC", "TG":"TG", "GT":"CA", "GC":"CG", "GA":"CT", "AT":"TA", "AG":"TC", "AC":"TG"} # for snp occurence substitution
sample_seq = []
line = dfh.readline().strip('\n').split('\t')
if line[8] in exonic_db.iterkeys():
	dict_sm[line[0]] = '1'
	dict_dsm[line[0]] += '/' + exonic_db[line[8]]

if len(line[11]) == 1 and len(line[12]) == 1:
	sig = sigsub[line[11] + line[12]] # if "TA" ==> "TA", "AT" ==> "TA"
	trinuc = comp(line[65][9:12].upper())
	detsig = trinuc + '_' + sig # detailed signature ex) ACT_CA
	sigdb[detsig] += 1
	titvdb[sig] += 1
sample = line[15]
sample_seq.append(sample)
tmbc = 1 # Mutation count for each sample for TMB (rfh3)
line = dfh.readline().strip('\n').split('\t')

#### Loop Start ####

while line != ['']:
	if line[15] == sample:
		tmbc += 1
		if line[8] in exonic_db.iterkeys():
			dict_sm[line[0]] = '1'
			dict_dsm[line[0]] += '/' + exonic_db[line[8]]
			if len(line[11]) == 1 and len(line[12]) == 1 and line[9] == "SNP":
			        sig = sigsub[line[11] + line[12]] # if "TA" ==> "TA", "AT" ==> "TA"
			        trinuc = comp(line[65][9:12].upper())
			        detsig = trinuc + '_' + sig # detailed signature ex) ACT_CA
			        sigdb[detsig] += 1
			        titvdb[sig] += 1
		else:
			if len(line[11]) == 1 and len(line[12]) == 1 and line[9] == "SNP":
                                sig = sigsub[line[11] + line[12]] # if "TA" ==> "TA", "AT" ==> "TA"
				trinuc = comp(line[65][9:12].upper())
				detsig = trinuc + '_' + sig # detailed signature ex) ACT_CA
                                sigdb[detsig] += 1
                                titvdb[sig] += 1
	else:
		tmb[sample] = tmbc
		unzipped_mutinfo.append(od_to_tuple(dict_sm))
		unzipped_dmutinfo.append(od_to_tuple(dict_dsm))
		rfh4.write(sample + '\t' + '\t'.join(map(lambda x: str(x), map(lambda x: x / float(sum(sigdb.itervalues())), sigdb.itervalues()))) + '\n')
		rfh5.write(sample + '\t' + '\t'.join(map(lambda x: str(x), map(lambda x: x / float(sum(titvdb.itervalues())), titvdb.itervalues()))) + '\n')
		rfh4.flush()
		rfh5.flush()
		dict_sm = OrderedDict(map(lambda x: (x, '0'), map(lambda y: y[0], sorted_t_mutc))) # dictionary of mutation presense of each sample
		dict_dsm = OrderedDict(map(lambda x: (x, '0'), map(lambda y: y[0], sorted_t_mutc))) # dictionary of "detailed" mutation presence of each sample
		sigdb = OrderedDict([("ACA_CA",0), ("ACC_CA",0), ("ACG_CA",0), ("ACT_CA",0), ("CCA_CA",0), ("CCC_CA",0), ("CCG_CA",0), ("CCT_CA",0), ("GCA_CA",0), ("GCC_CA",0), ("GCG_CA",0), ("GCT_CA",0), ("TCA_CA",0), ("TCC_CA",0), ("TCG_CA",0), ("TCT_CA",0), ("ACA_CG",0), ("ACC_CG",0), ("ACG_CG",0), ("ACT_CG",0), ("CCA_CG",0), ("CCC_CG",0), ("CCG_CG",0), ("CCT_CG",0), ("GCA_CG",0), ("GCC_CG",0), ("GCG_CG",0), ("GCT_CG",0), ("TCA_CG",0), ("TCC_CG",0), ("TCG_CG",0), ("TCT_CG",0), ("ACA_CT",0), ("ACC_CT",0), ("ACG_CT",0), ("ACT_CT",0), ("CCA_CT",0), ("CCC_CT",0), ("CCG_CT",0), ("CCT_CT",0), ("GCA_CT",0), ("GCC_CT",0), ("GCG_CT",0), ("GCT_CT",0), ("TCA_CT",0), ("TCC_CT",0), ("TCG_CT",0), ("TCT_CT",0), ("ATA_TA",0), ("ATC_TA",0), ("ATG_TA",0), ("ATT_TA",0), ("CTA_TA",0), ("CTC_TA",0), ("CTG_TA",0), ("CTT_TA",0), ("GTA_TA",0), ("GTC_TA",0), ("GTG_TA",0), ("GTT_TA",0), ("TTA_TA",0), ("TTC_TA",0), ("TTG_TA",0), ("TTT_TA",0), ("ATA_TC",0), ("ATC_TC",0), ("ATG_TC",0), ("ATT_TC",0), ("CTA_TC",0), ("CTC_TC",0), ("CTG_TC",0), ("CTT_TC",0), ("GTA_TC",0), ("GTC_TC",0), ("GTG_TC",0), ("GTT_TC",0), ("TTA_TC",0), ("TTC_TC",0), ("TTG_TC",0), ("TTT_TC",0), ("ATA_TG",0), ("ATC_TG",0), ("ATG_TG",0), ("ATT_TG",0), ("CTA_TG",0), ("CTC_TG",0), ("CTG_TG",0), ("CTT_TG",0), ("GTA_TG",0), ("GTC_TG",0), ("GTG_TG",0), ("GTT_TG",0), ("TTA_TG",0), ("TTC_TG",0), ("TTG_TG",0), ("TTT_TG", 0)])
		titvdb = OrderedDict([("CT",0),("TC",0),("TG",0), ("TA",0), ("CG",0), ("CA",0)])
		if line[8] in exonic_db.iterkeys():
			dict_sm[line[0]] = '1'
			dict_dsm[line[0]] += '/' + exonic_db[line[8]]
			if len(line[11]) == 1 and len(line[12]) == 1 and line[9] == "SNP":
                                sig = sigsub[line[11] + line[12]] # if "TA" ==> "TA", "AT" ==> "TA"
                                trinuc = comp(line[65][9:12].upper())
                                detsig = trinuc + '_' + sig # detailed signature ex) ACT_CA
                                sigdb[detsig] += 1
                                titvdb[sig] += 1
		else:	
			if len(line[11]) == 1 and len(line[12]) == 1 and line[9] == "SNP":
				sig = sigsub[line[11] + line[12]] # if "TA" ==> "TA", "AT" ==> "TA"
				trinuc = comp(line[65][9:12].upper())
				detsig = trinuc + '_' + sig # detailed signature ex) ACT_CA
				sigdb[detsig] += 1
				titvdb[sig] += 1
		sample = line[15]
		sample_seq.append(sample)
		tmbc = 1
	line = dfh.readline().strip('\n').split('\t')
tmb[sample] = tmbc
unzipped_mutinfo.append(od_to_tuple(dict_sm))
unzipped_dmutinfo.append(od_to_tuple(dict_dsm))
rfh4.write(sample + '\t' + '\t'.join(map(lambda x: str(x), map(lambda x: x / float(sum(sigdb.itervalues())), sigdb.itervalues()))) + '\n')
rfh5.write(sample + '\t' + '\t'.join(map(lambda x: str(x), map(lambda x: x / float(sum(titvdb.itervalues())), titvdb.itervalues()))) + '\n')
final = zip(*unzipped_mutinfo) # final matrix
final_d = zip(*unzipped_dmutinfo) # final detailed matrix

# Write Final Output
rfh1.write('Gene\t' + '\t'.join(sample_seq) + '\n')
for i in range(len(dict_stm)):
	rfh1.write(tuple(dict_stm)[i] + '\t' + '\t'.join(list(final[i])) + '\n')
	rfh1.flush()
rfh2.write('Gene\t' + '\t'.join(sample_seq) + '\n')
for i in range(len(dict_stm)):
        rfh2.write(tuple(dict_stm)[i] + '\t' + '\t'.join(list(final_d[i])) + '\n')
        rfh2.flush()
rfh3.write('sample\tTMB(Total)\n')
for i in sample_seq:
	rfh3.write(i + '\t' + str(tmb[i]) + '\n')
	rfh3.flush()
dfh.close()
rfh1.close()
rfh2.close()
rfh3.close()
dbf = open("CosmicCodingMutationCensus.txt", 'r')
cosmicdb = []
for i in dbf:
	cosmicdb.append(i.strip('\n').split()[0])
dfh1 = open(sys.argv[1][:-7] + '_matrix.txt', 'r') # 0: no mut, 1: yes mut
dfh2 = open(sys.argv[1][:-7] + '_matrix_detail.txt', 'r') # detailed mutational classfication, TMB (filter criteria needed), mutational spectrum...
rfh1 = open('1.txt', 'w')
rfh2 = open('2.txt', 'w')
rfh1.write(dfh1.readline())
rfh2.write(dfh2.readline())
for i in dfh1:
	line = i.strip('\n').split('\t')
	if line[0] in cosmicdb:
		rfh1.write('\t'.join(line) + '\n')
		rfh1.flush()
for i in dfh2:
	line = i.strip('\n').split('\t')
	if line[0] in cosmicdb:
		rfh2.write('\t'.join(line) + '\n')
		rfh2.flush()
dfh1.close()
dfh2.close()
rfh1.close()
rfh2.close()
dfh1 = open("1.txt", 'r')
rfh1 = open("3.txt", 'w')
fuck = list(i.strip('\n').split('\t') for i in dfh1)
for i in zip(*fuck):
	rfh1.write('\t'.join(i) + '\n')
	rfh1.flush()
dfh1.close()
rfh1.close()
dfh2 = open("2.txt", 'r')
rfh2 = open("4.txt", 'w')
shit = list(i.strip('\n').split('\t') for i in dfh2)
for i in zip(*shit):
        rfh2.write('\t'.join(i) + '\n')
        rfh2.flush()
dfh2.close()
rfh2.close()
call('/bin/rm -rf 1.txt', shell=True)
call('/bin/rm -rf 2.txt', shell=True)
call('/bin/mv -f 3.txt ' + sys.argv[1][:-7] + '_matrix.txt', shell=True)
call('/bin/mv -f 4.txt ' + sys.argv[1][:-7] + '_matrix_detail.txt', shell=True)

dfh = open(sys.argv[1][:-7] + '_matrix.txt', 'r')
tl = len(dfh.readline().strip('\n').split('\t')) # table length for linux sorting
sort_command = map(lambda y: str(y) + ',' + str(y) + 'r', range(2, tl + 1))
dfh.close()
call('/bin/sort -k ' + ' -k '.join(sort_command) + ' ' + sys.argv[1][:-7] + '_matrix.txt > fuck.txt', shell=True)

dfh = open("fuck.txt", 'r')
rfh = open(sys.argv[1][:-7] + '_matrix_finale.txt', 'w')
fuck = list(i.strip('\n').split('\t') for i in dfh)
for i in zip(*fuck):
	rfh.write('\t'.join(i) + '\n')
	rfh.flush()
dfh.close()
rfh.close()
call('/bin/rm -rf fuck.txt', shell=True)

dfh = open(sys.argv[1][:-7] + '_matrix_detail.txt', 'r')
shit = dfh.readline().strip().split()
det_db = {} # detailed DB
for i in range(len(shit)-1):
	for j in dfh:
		line = j.strip().split()
		try:
			det_db[line[0]] += '|' + line[i+1]
		except KeyError:
			det_db[line[0]] = line[i+1]
	dfh.seek(0)
	dfh.readline()
dfh.close()
dfh = open(sys.argv[1][:-7] + '_matrix_finale.txt', 'r')
rfh = open(sys.argv[1][:-7] + '_matrix_finale_detail.txt', 'w')
header = dfh.readline().strip().split()[1:]
rfh.write('Gene\t' + '\t'.join(header) + '\n')
c = 0
for i in dfh:
	line = i.strip().split('\t')
	rfh.write(line[0] + '\t' + '\t'.join(map(lambda x: det_db[header[x]].split('|')[c], range(len(header)))) + '\n')
	c += 1
dfh.close()
rfh.close()
call('/bin/rm -rf ' + sys.argv[1][:-7] + '_matrix_finale.txt', shell=True)
