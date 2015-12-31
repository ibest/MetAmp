import sys
import collections
from Bio import SeqIO


#../usearch/usearch7.0.1090_i86osx32 -acceptall -allpairs_global gold1500.fa -uc gold1500.dmx

inhandle = open(sys.argv[-1],'rU') #gold1500.dmx
outhandle = open("duplicated.txt",'w')
dup_ids = []
for line in inhandle.readlines() :
	splitted = line.split('\t')
	if splitted[0] != 'N' :
		if float(splitted[3]) > 99.00 :
			dup_id = splitted[9].split('/')[1]
			#print dup_id
			dup_ids.append(dup_id.strip())
			outhandle.write(dup_id+'\n')
		
outhandle.close()
inhandle.close()
print len(dup_ids)
print [x for x, y in collections.Counter(dup_ids).items() if y > 2]

#print dup_ids
inhandle = open("gold1500.fa",'rU')
outhandle = open("gold1500_no_dup.fa",'w')
i=0
for record in SeqIO.parse(inhandle,"fasta") :
	if record.id.split('/')[1] not in dup_ids :
		i+=1
		SeqIO.write(record,outhandle,'fasta')

outhandle.close()
print i, len(dup_ids)
inhandle.close()

