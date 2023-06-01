import os, sys
from Bio import SeqIO
import re

seqs = list(SeqIO.parse(sys.argv[1], 'fasta'))

unique = set()
temp = []
for seq in seqs:
	if seq.id not in unique:
		unique.add(seq.id)
		temp.append(seq)

seqs = temp

seqs = [x for x in seqs if x.seq.count("N") / len(x) < 0.1]

newseqs = []
for x in seqs:
	search1 = re.search("([^\s]+)[^\(]*\(([^\)]+\))", x.description)
	search2 = re.search('(H\dN\d)', x.description)
	if search1 and search2:
		x.id = "|".join(search1.groups()).replace(" ","_") + "|A|" + search2.group(1) 
		x.name = x.id
		x.description = ''
		newseqs.append(x)

newseqs = [x for x in newseqs if "|A/" in x.id]


SeqIO.write(newseqs, sys.argv[2], 'fasta')
