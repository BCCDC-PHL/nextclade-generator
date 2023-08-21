#%%
import os, sys
from Bio import SeqIO
import argparse
import multiprocessing as mp
from tools import pairwise_alignment, flex_translate, init_aligner

os.chdir(os.path.dirname(__file__))

#%%
def init_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s",'--sequences', required=True, help='Path to a FASTA file containing all sequences that will be participating in the dataset')
	parser.add_argument("-r",'--reference', default='../assets/h5_ref.fasta', help='Path to a FASTA file containing all sequences that will be participating in the dataset')
	parser.add_argument("-m",'--metadata', required=True, help='Path to a metadata describing all samples participating in the dataset')
	parser.add_argument("-t",'--tree', required=True, help='Path to a metadata describing all samples participating in the dataset')
	parser.add_argument("-n",'--threads', default=1, type=int, help='Number of threads to use when computing mutations.')

def get_nucl_mutations(aligner, ref, qry):
	ref_aln, qry_aln = pairwise_alignment(aligner, ref, qry)
	
	mut_list = []
	for n, (char1, char2) in enumerate(zip(ref_aln, qry_aln),1):
		if char1 != char2 and char1 not in ['*','-','?'] and char2 not in ['*','-','?']:
			mut_list.append((int(n), char2))
	return mut_list

def transform(pos):
	pos = int(pos)
	if pos > 353: 
		return "HA2", pos - 353
	elif pos > 23:
		return "HA1", pos - 23
	elif pos > 7:
		return "SigPep", pos - 7
	else:
		return "ERROR", -1

def get_aa_mutations(aligner, aa_ref, nt_qry):
	aa_qry, _, _ = flex_translate(nt_qry)

	ref_aln, qry_aln = pairwise_alignment(aligner, aa_ref, aa_qry)

	mut_list = []
	for n, (char1, char2) in enumerate(zip(ref_aln, qry_aln),1):
		if char1 != char2 and char1 not in ['*','-','?'] and char2 not in ['*','-','?']:
			gene, position = transform(n)
			if gene != "ERROR" and gene != "SigPep":
				mut_list.append((gene, int(position), char2))

	return mut_list

#%%
def get_all_mutations(ref_path, qry_path, nthreads):

	nuc_ref = next(SeqIO.parse(ref_path,'fasta'))
	

	aa_ref = nuc_ref.translate()
	nuc_ref = str(nuc_ref.seq)

	seqs = list(SeqIO.parse(qry_path, 'fasta'))
	aligner = init_aligner()

	global nucl_align_wrap    		# necessary to run function in parallel

	def nucl_align_wrap(x):
		return get_nucl_mutations(aligner, nuc_ref, str(x.seq))
	

	global aa_align_wrap 			# necessary to run function in parallel

	def aa_align_wrap(x):
		return get_aa_mutations(aligner, aa_ref, str(x.seq))


	print("Starting parallelized.")
	with mp.Pool(nthreads) as pool:
		ids = [x.id for x in seqs]
		nucl_dict = list(pool.map(nucl_align_wrap, seqs))
		nucl_dict = dict(zip(ids,nucl_dict))
		aa_dict = list(pool.map(aa_align_wrap, seqs))
		aa_dict = dict(zip(ids,aa_dict))

	print("Finishing parallelized.")

	return nucl_dict, aa_dict

if __name__ == "__main__":

	parser = init_parser()
	args = parser.parse_args()

	nucl_dict, aa_dict = get_all_mutations(args.reference, args.sequences, args.nthreads)