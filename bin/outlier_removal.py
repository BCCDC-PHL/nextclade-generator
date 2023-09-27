import os, sys 
import pandas as pd 
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


#%%
def make_consensus(alignment):
	centroid_seq = []

	for i in range(alignment.get_alignment_length()):
		# Get the most common character/residue in each column of the alignment
		column = alignment[:, i]
		most_common_residue = max(set(column), key=column.count)
		centroid_seq.append(most_common_residue)

	centroid_seq = "".join(centroid_seq)

	centroid = SeqRecord(
		id='CONSENSUS',
		name='CONSENSUS',
		description='CONSENSUS',
		seq=Seq(centroid_seq)
	)

	return centroid

def compute_distances(centroid, alignment):

	calculator = DistanceCalculator('identity')
	alignment.extend([centroid])

	distances = calculator.get_distance(alignment)

	return distances

#%%
# Load sequences with SeqIO
base_path = '/home/john.palmer/work/flu/testing/nc_gen/'
sequences = list(SeqIO.parse(base_path + "master.fa", "fasta"))
#%%
# Perform Multiple Sequence Alignment using MUSCLE
alignment = AlignIO.read(base_path + "master_align.fa", "fasta")

centroid = make_consensus(alignment)

#%% 
dist_mat = compute_distances(centroid, alignment)

dist_df = pd.DataFrame(list(dist_mat), index=dist_mat.names, columns=dist_mat.names)

distances = dist_df.loc['CONSENSUS']
distances = distances.iloc[0:-1]
# %%
# compute the centroid of each cluster, where each cluster is a clade 

# calculate the pairwise distances of each sequence to their centroid

# those clusters with a STDEV higher than a certain extent will be subject to pruning
# since this is indicative of having weird outliers 

# sequences with a STDEV with magnitude > 2.0 will be removed from the dataset