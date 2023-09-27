from Bio import SeqIO, Entrez
import os, sys
import pandas as pd 
import numpy as np


# environment variables need to be set with NCBI credentials, both email and API key (from creating an account)
Entrez.email = os.environ['NCBI_EMAIL']
Entrez.api_key = os.environ['NCBI_API_KEY']


handle = Entrez.esearch(db="nucleotide", retmax=10, term="(A/HongKong/485/1997) AND (HA)", idtype="acc")
record = Entrez.read(handle)
handle.close()

print(record['Count'])
print(record['IdList'])


def sample_clades(filepath):
	
	clade_dict = {}
	with open(filepath, 'r') as infile:
		infile.readline()
		for n, line in enumerate(infile.readlines()):
			if n % 100 == 0:
				print(n)
			name, clade = line.split()

			if clade not in clade_dict:
				clade_dict[clade] = []
			clade_dict[clade].append(name)
	if "?" in clade_dict:
		del clade_dict['?']

	for clade in clade_dict.keys():
		if len(clade_dict[clade]) > 100 :
			clade_dict[clade] = np.random.choice(clade_dict[clade], 100).tolist()
		
	final_accnos = [y for x in clade_dict for y in clade_dict[x]]
	return final_accnos

def fetch_ids(accno_list):

	accnos = []
	num_added = []
	outfile = open('h5n1_clades.fasta', 'w')

	for n, accno in enumerate(accno_list):
		if n % 100 == 0 :
			print(n)
	
		handle = Entrez.esearch(db="nucleotide", retmax=5, term=f"({accno}) AND (HA)", idtype="acc")
		record = Entrez.read(handle)
		handle.close()
		if int(record['Count']) > 0:
			accnos += record['IdList']
			num_added = record['Count']

			query = Entrez.efetch(db='nucleotide', id=record['IdList'][0], rettype='fasta', retmode='text')

			for i in query:
				outfile.write(i)

	outfile.close()

	return True

def fetch_entries(accno_list, db_name, return_type, outpath):

	accnos = pd.Series(accno_list)

	accnos = accnos.drop_duplicates()
	
	print(accnos[0:10])

	query = Entrez.efetch(db='nucleotide', id=accnos, rettype='fasta', retmode='text')

	with open(outpath, 'w') as outfile:
		for i in query:
			outfile.write(i)


if __name__ == '__main__':
	clade_list =  sample_clades(sys.argv[1])
	results = fetch_ids(clade_list)
	
