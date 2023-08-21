import os 
import sys
from collections import defaultdict
import re
import numpy as np
from Bio import Align
from Bio.SeqRecord import SeqRecord


#%%
def parse_fasta(filepath):	
	seqs = {}
	with open(filepath, 'r') as handle:
		for line in handle.readlines():
			if line[0] == '>':
				header = line.strip().lstrip('>')
				seqs[header] = ''
			else:
				seqs[header] += line.strip()
	return seqs	

#%%
def write_fasta(seqs, outpath):
	try: 
		with open(outpath, 'w') as outfile:

			for header, seq in seqs.items():
				outfile.write(">" + header + '\n')
				outfile.write(seq + '\n')
	except Exception as e:
		print(str(e))
		return False
	return True

complement_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 
                    'W':'S', 'R':'Y', 'K':'M', 'Y':'R', 'S':'W', 'M':'K',
                    'B':'V', 'D':'H', 'H':'D', 'V':'B',
                    '*':'*', 'N':'N', '-':'-'}

def reverse_and_complement(seq):
    rseq = seq[::-1]
    rcseq = ''
    for i in rseq:  # reverse order
        rcseq += complement_dict[i]
    return rcseq

codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 
'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
'---':'-', 'XXX':'?'}

mixture_regex = re.compile('[WRKYSMBDHVN-]')

mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT', 'S':'CG', 
'M':'AC', 'V':'AGC', 'H':'ATC', 'D':'ATG', 
'B':'TGC', 'N':'ATGC', '-':'ATGC'}

ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.items())


def translate_nuc(seq, offset, resolve=False, return_list=False):
	"""
	Translate nucleotide sequence into amino acid sequence.
		offset by X shifts sequence to the right by X bases
	Synonymous nucleotide mixtures are resolved to the corresponding residue.
	Nonsynonymous nucleotide mixtures are encoded with '?' 
	"""
	
	seq = '-'*offset + seq
	
	aa_list = []
	aa_seq = ''	# use to align against reference, for resolving indels
	
	# loop over codon sites in nucleotide sequence
	for codon_site in range(0, len(seq), 3):
		codon = seq[codon_site:codon_site+3]
		
		if len(codon) < 3:
			break
		
		# note that we're willing to handle a single missing nucleotide as an ambiguity
		if codon.count('-') > 1 or '?' in codon:
			if codon == '---':	# don't bother to translate incomplete codons
				aa_seq += '-'
				aa_list.append(['-'])
			else:
				aa_seq += '?'
				aa_list.append(['?'])
			continue
		
		# look for nucleotide mixtures in codon, resolve to alternative codons if found
		num_mixtures = len(mixture_regex.findall(codon))
		
		if num_mixtures == 0:
			aa = codon_dict[codon]
			aa_seq += aa
			aa_list.append([aa])
			
		elif num_mixtures == 1:
			resolved_AAs = []
			for pos in range(3):
				if codon[pos] in mixture_dict.keys():
					for r in mixture_dict[codon[pos]]:
						rcodon = codon[0:pos] + r + codon[(pos+1):]
						if codon_dict[rcodon] not in resolved_AAs:
							resolved_AAs.append(codon_dict[rcodon])
							
			aa_list.append(resolved_AAs)
			
			if len(resolved_AAs) > 1:
				if resolve:
					# for purposes of aligning AA sequences
					# it is better to have one of the resolutions
					# than a completely ambiguous '?'
					aa_seq += resolved_AAs[0]
				else:
					aa_seq += '?'
			else:
				aa_seq += resolved_AAs[0]
				
		else:
			aa_seq += '?'
			aa_list.append(['?'])
			
	if return_list:
		return aa_list

	return aa_seq

def flex_translate(nt_seq, debug=False):
	"""
	Function to find the most appropriate amino acid translation by testing all three reading frames.
	Designed to receive Biopython SeqRecord objects.
	Returns three outputs as a tuple:
		- A SeqRecord object containing the best amino acid translation 
		- An integer indicating the best reading frame that was used for translation (0, 1, or 2)
		- An integer indicating the number of stop codons in this best translation candidate
	"""
	min_count = np.Inf
	best_seq = None
	best_frame = -1

	def pad(seq):
		"""
		Supporting function to prevent Biopython translate() function from complaining about non-multiples-of-3
		"""
		mod = len(seq) % 3
		if mod != 0:
			toadd = 3 - mod
			seq += "N" * toadd
		return seq

	def translate_function(nt_seq):
		"""
		Allows for type agnostic amino acid translation. 
		Works for both strings and SeqRecord objects. Benefit is that SeqRecord objects are maintained as SeqRecords.
		"""
		if isinstance(nt_seq, SeqRecord):
			return nt_seq.translate()
		else:
			return translate_nuc(nt_seq, 0)

	for i in range(3):
		nt_tmp = pad(nt_seq[i:])
		aa_seq = translate_function(nt_tmp)
		aa_count = aa_seq.count("*")
		
		if debug: 
			print(aa_seq)
	
		if aa_count < min_count:
			min_count = aa_count
			best_seq = aa_seq
			best_frame = i

	if isinstance(best_seq, SeqRecord):
		best_seq.id = nt_seq.id
		best_seq.name = nt_seq.name
		best_seq.description = nt_seq.description
	
	return best_seq, best_frame, min_count

def init_aligner(mode='global', open_gap=-1.0, x_gap=-0.1):
	aligner = Align.PairwiseAligner()
	aligner.mode = mode
	aligner.open_gap_score = open_gap
	aligner.extend_gap_score = x_gap
	aligner.target_end_gap_score = 0.0
	aligner.query_end_gap_score = 0.0
	return aligner

def get_boundaries(str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')
    # return a tuple giving indices of subsequence without gap prefix and suffix
    res = [0,len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])
    if right:
        res[1] = len(str) - len(right[0])
    return res

# performs a pairwise alignment between two Biopython SeqRecord objects
def pairwise_alignment(aligner, ref, qry):
	"""
	Arguments
	- `aligner`: an object of the Biopython alignment tool to be used for alignment
	- `ref`: a Bio.SeqRecord object representing the reference sequence
	- `qry`: a Bio.SeqRecord object representing the query sequence

	Returns
	- `ref_aln`: a string representing the aligned reference sequence
	- `qry_aln`: a string representing the aligned query sequence

	The function performs a pairwise alignment between the reference and query sequences using the provided `aligner` object from the Biopython alignment tool. 
	The resulting alignment is then cut down to only the region of interest that overlaps with the reference sequence, as determined by the `get_boundaries()` function. 
	The aligned reference and query sequences are returned as strings.
	"""
	if isinstance(ref, SeqRecord):
		ref = str(ref.seq)
	if isinstance(qry, SeqRecord):
		ref = str(ref.seq)
	ref_aln, qry_aln = next(aligner.align(ref, qry))
	# cut down the alignment to only the region of interest (that which overlaps with reference)
	start, stop = get_boundaries(ref_aln)
	ref_aln = ref_aln[start:stop]
	qry_aln = qry_aln[start:stop]
	return ref_aln, qry_aln