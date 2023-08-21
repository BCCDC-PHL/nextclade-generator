import os, sys
from Bio import SeqIO
import re

seqs = list(SeqIO.parse(sys.argv[1], 'fasta'))

def filter_unique(seq_list):
	unique = set()
	filtered = []
	for seq in seq_list:
		if seq.id not in unique:
			unique.add(seq.id)
			filtered.append(seq)

	return filtered


def filter_N(seq_list):
	return [x for x in seq_list if x.seq.count("N") / len(x) < 0.1]

def rename_headers(seq_list):
	newseqs = []
	for x in seq_list:
		re_fields = re.search("([^\s]+)[^\(]*\(([^\)]+\))", x.description)
		re_subtype = re.search('(H\dN\d)', x.description)

		if re_fields and re_subtype:
			x.id = "|".join(re_fields.groups()).replace(" ","_") + "|A|" + re_subtype.group(1) 
			x.name = x.id
			x.description = ''
			newseqs.append(x)
	
	return [x for x in newseqs if "|A/" in x.id]


seqs = filter_unique(seqs)
seqs = filter_N(seqs)
seqs = rename_headers(seqs)


SeqIO.write(seqs, sys.argv[2], 'fasta')
