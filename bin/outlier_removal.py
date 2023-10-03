#%%
import os, sys 
import pandas as pd 
import numpy as np 
import json
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import run, PIPE
import multiprocessing as mp 
from collections import defaultdict
from tools import get_mutations, get_boundaries
from sklearn.ensemble import IsolationForest
#%%
def count_column(column):
	return max(set(column), key=column.count)


def make_consensus(alignment, clade_name):

	centroid_seq = []

	alignment_cols = [alignment[:,i] for i in range(alignment.get_alignment_length())]

	if len(alignment) > 300:
		threads = len(alignment) // 100
		threads = min(threads, 32)
		print(f"Making consensus in parallel with {threads} threads...")


		with mp.Pool(processes=threads) as pool:
			centroid_seq = list(pool.map(count_column, alignment_cols, chunksize=100))
	
	else:
		print('Making consensus...')
		centroid_seq = [max(set(column), key=column.count) for column in alignment_cols]

	# for i in range(alignment.get_alignment_length()):
	# 	# Get the most common character/residue in each column of the alignment
	# 	column = alignment[:, i]
	# 	most_common_residue = max(set(column), key=column.count)
	# 	centroid_seq.append(most_common_residue)

	centroid_seq = "".join(centroid_seq)

	name = 'CONSENSUS|' + clade_name

	centroid = SeqRecord(
		id=name,
		name=name,
		description=name,
		seq=Seq(centroid_seq)
	)

	return centroid

#%%
def compute_distances(centroid, alignment):

	calculator = DistanceCalculator('identity')
	alignment.extend([centroid])

	distances = calculator.get_distance(alignment)

	return distances

def make_consensus_seqs(seq_path, meta_path, meta_delim='\t'):
	
	os.makedirs('alignments', exist_ok=True)

	if os.path.isfile('consensus_seqs.fasta'):
		return list(SeqIO.parse('consensus_seqs.fasta', 'fasta'))

	metadata = pd.read_csv(meta_path, sep=meta_delim)
	clade_dict = metadata.groupby("clade")['name'].apply(set).to_dict()

	consensus_seqs = []

	for clade, ids in clade_dict.items():

		if clade.endswith("-like"):
			continue

		print(clade)

		clade_name = clade.replace(".","_")

		tmp_seq_path = "alignments/.temp.fasta"
		out_aln_path = 'alignments/' + clade_name + '.aln.fasta'

		if not os.path.isfile(out_aln_path):
			seqs = [seq for seq in SeqIO.parse(seq_path, "fasta") if seq.id in ids]

			SeqIO.write(seqs, tmp_seq_path, 'fasta')

			threads = len(seqs) // 100
			threads = max(threads, 4)
			threads = min(threads, 32)

			print(f"Starting MAFFT with {threads} threads...")

			result = run(f"mafft --auto --thread {threads} {tmp_seq_path} > {out_aln_path}", shell=True, capture_output=True)

			if result.returncode != 0:
				print("ERROR: Mafft alignment process failed. ")
				print(result.stdout.decode('utf-8'))
				print(result.stderr.decode('utf-8'))
				exit(1)

			print("Completed MAFFT successfully.")

		alignment = AlignIO.read(out_aln_path, "fasta")

		consensus = make_consensus(alignment, clade_name)

		consensus_seqs += (consensus,)
	
	SeqIO.write(consensus_seqs, 'consensus_seqs.fasta', 'fasta')

	return consensus_seqs

#%%
def outlier_pairwise(consensus_seqs):

	distance_dict = {}

	clades = [seq.id.split("|")[-1] for seq in consensus_seqs]

	for clade_name, consensus_seq in zip(clades, consensus_seqs):
		print(clade_name)

		alignment = AlignIO.read('alignments/' + clade_name + '.aln.fasta', "fasta")

		dist_mat = compute_distances(consensus_seq, alignment)

		dist_df = pd.DataFrame(list(dist_mat), index=dist_mat.names, columns=dist_mat.names)

		distances = dist_df.iloc[-1,0:-1]


		upper_qrt, lower_qrt = np.quantile(distances, (0.75, 0.25)) 
		iqr = upper_qrt - lower_qrt
		stdev = np.std(distances)

		upper_bound =  upper_qrt + 2.5 * iqr

		outlier_mask = (distances > upper_bound)
		outliers = distances.loc[outlier_mask]
		non_outliers = distances.loc[~outlier_mask]

		print(f"IQR : {iqr}")
		print(f"SD : {stdev}")
		print(f"Upper Bound : {upper_bound}")
		print(f"Outliers : {outliers}")
		print(f"Non-outliers : {non_outliers.head()}")

		distance_dict[clade_name] = distances


	return distance_dict


#%%

def outlier_isolation_forest(consensus_seqs):


	clades = [seq.id.split("|")[-1] for seq in consensus_seqs]

	for clade_name, consensus_seq in zip(clades, consensus_seqs):
		print(clade_name)

		alignment = AlignIO.read('alignments/' + clade_name + '.aln.fasta', "fasta")

		start, stop = get_boundaries(consensus_seq)

		align_cut = alignment[:, start:stop]

		consensus_seq = consensus_seq[start:stop]

		mutation_dict = {} 

		for seq in align_cut:
			mutation_dict[seq.id] = get_mutations(consensus_seq, seq)

		all_mutations = []

		for name in mutation_dict:
			for mut in mutation_dict[name]:
				if len(mut[0])==1 and len(mut[2])==1 and mut[2]!='-':
					all_mutations.append((name, mut[0]+str(mut[1])+mut[2]))

		mutations_df = pd.DataFrame(all_mutations,  columns=['name','mutation'])
		mutations_df['present'] = 1
		mutations_df = pd.pivot_table(mutations_df, values='present', columns='mutation', index='name').fillna(0)



		break
		
	return mutation_dict


def valid_ascent(previous, current):
	if current.startswith(previous):
		return True



#%%
base_path = '/home/john.palmer/work/flu/testing/nc_gen/'
os.chdir(base_path)

consensus_seqs = make_consensus_seqs(base_path + "master.fasta", base_path + "master_meta.tsv")


#%%
distance_dict = outlier_pairwise(consensus_seqs)

#%%
outlier_isolation_forest(consensus_seqs)
# %%
# compute the centroid of each cluster, where each cluster is a clade 

# calculate the pairwise distances of each sequence to their centroid

# those clusters with a STDEV higher than a certain extent will be subject to pruning
# since this is indicative of having weird outliers 

# sequences with a STDEV with magnitude > 2.0 will be removed from the dataset

#%%

# def find_max_mono(node, mono_clade_dict):
# 	if 'children' not in node:
# 		node_clade = node['node_attrs']['clade_membership']['value'] if 'clade_membership' in node['node_attrs'] else 'NA'
# 		return [node['name']], [node_clade]

# 	child_names = []
# 	child_clades = []
# 	for child in node.get('children', []):
# 		name_list, clade_list = find_max_mono(child, mono_clade_dict)
# 		child_names += name_list
# 		child_clades += clade_list


# 	if len(set(child_clades)) == 1:
# 		clade = child_clades[0]

# 		if clade not in mono_clade_dict or len(child_names) > len(mono_clade_dict[clade]):
# 			mono_clade_dict[clade] = child_names

# 	return child_names, child_clades


# #%% augur tree filtering and removal
# base_path = '/home/john.palmer/git/nextclade-generator/nextclade/'

# #treepath = '/home/john.palmer/work/flu/ref/augur_h5n1/ha_dataset/master2/auspice/tree.json'
# treepath = base_path + 'auspice/tree_tip_label.json'
# with open(treepath, 'r') as handle:
# 	data = json.load(handle)
# tree = data['tree']

# %%
