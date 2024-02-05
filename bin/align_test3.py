#%%
import os, sys
from Bio import Phylo
import pandas as pd
import numpy as np 
from collections import defaultdict, Counter
from functools import partial
#import augur_utils
import importlib
from subprocess import run
import json
from compute_mutations import get_all_mutations

os.chdir(os.path.dirname(__file__))
pd.set_option('display.max_rows',999)
os.makedirs('../nextclade/clades',exist_ok=True)

#%%
base_path = '/home/john.palmer/git/nextclade-generator/nextclade2/'
input_seq_path = base_path + 'inputs/master-final-nodup.fasta'
input_metadata_path = base_path + "inputs/master-meta-final-nodup.tsv"
input_tree = base_path + 'refine/h5_refined.nwk'
aa_json_path = base_path + 'translate/h5_aa_node_data.json'
ancestral_json_path = base_path + 'ancestral/h5_anc_node_data.jsonn'
nthreads = 6
nextstrain_env = '/home/john.palmer/miniconda/envs/nextstrain'

#%%
def tune_thresholds(thres_dict, clade_dict, nucl_dict, aa_dict):
	print("Computing nucleotide and amino acid mutations...")
	final_mutations = {}

	for clade in clade_dict:
		
		N = len(clade_dict[clade])
		aa_thres = N * 0.95 
		nucl_thres = N * 0.95 

		if clade in thres_dict:
			aa_thres, nucl_thres = thres_dict[clade] * N, thres_dict[clade] * N, 

		if clade in "8 0 ".split():
			continue

		master_aa_list = [aa for id in clade_dict[clade] for aa in aa_dict[id]]
		master_nuc_list = [nuc for id in clade_dict[clade] for nuc in nucl_dict[id]]
		
		aa_counter = Counter(master_aa_list)
		nuc_counter = Counter(master_nuc_list)

		final_mutations[clade] = [x for x,y in aa_counter.items() if y > aa_thres]
		final_mutations[clade] += [('nuc',*x) for x,y in nuc_counter.items() if y > nucl_thres]

	
	final_mutations_list = [(clade,*x) for clade, vals in final_mutations.items() for x in vals]

	print("Finished computing mutations.")

	return final_mutations_list

def write_clades_tsv(final_mutations_list):
	print("Reformatting final mutations...")
	master_df = pd.DataFrame(final_mutations_list, columns=['clade','gene','site','alt'])
	master_df = master_df.sort_values(['clade','gene','site']).reset_index(drop=True)
	master_df['site'] = master_df['site'].astype(int)
	master_df.to_csv('../nextclade/clades/h5_clades.tsv', sep='\t', index=False)
	print("Final mutations written to h5_clades.tsv.")


def run_augur_clades():
	print("Starting augur clades...")
	result = run(f"cd ../nextclade; \
		{nextstrain_env}/bin/augur clades \
		--tree refine/h5_refined.nwk  \
		--mutations ancestral/h5_anc_node_data.json translate/h5_aa_node_data.json \
		--clades clades/h5_clades.tsv \
		--output-node-data clades/clades.json", shell=True)
	print("Finished augur clades.")

# --reference ../assets/h5n1_ref.fasta \
#%%
def compute_performance(meta):
	print("Starting to compute performance ...")

	jsondict = json.load(open('../nextclade/clades/clades.json'))
	predictions = [(id,jsondict['nodes'][id]['clade_membership']) for id in jsondict['nodes']]
	predictions = pd.DataFrame(predictions, columns=['name','pred'])
	predictions = predictions.loc[~predictions['name'].str.startswith("NODE")]
	predictions = predictions.merge(meta, on='name')

	N = predictions.shape[0]

	check = (predictions['pred'] == predictions['clade'])

	print(check.sum() / N)

	##%% calling a parent clade properly, but not the perfect match 
	partial = predictions[['pred', 'clade']].apply(lambda x: x['clade'].startswith(x['pred']), axis=1)
	print(partial.sum() / N)


	wrong = predictions.loc[~check]
	wrong = wrong.sort_values(['pred','clade'])
	print(wrong[['pred','clade']].value_counts())
	print(wrong['clade'].value_counts())
	print(wrong['pred'].value_counts())

	# update_parameters(predictions, thres_dict, best_thres, best_score, "clade", 'predictions')
	
	# track the maximum of each clade 
	# if the clade's performance has not improved over the past 10 iterations, freeze that clade's threshold in place 

	print("Finished performance computation.")
	


#%% runs a parallelized analysis that retrieves all nucleotide and amino acid mutations 
nucl_dict, aa_dict = get_all_mutations("../assets/h5n1_ref.fasta", input_seq_path, nthreads)

#%% find the intersection of common mutations 
meta = pd.read_csv(input_metadata_path,sep='\t')
meta = meta[['name','clade']]
clade_dict = meta.groupby('clade')['name'].apply(list).to_dict()

#%%
thres_dict = {clade : 0.95 for clade in clade_dict}

#%%
final_mutations_list = tune_thresholds(thres_dict, clade_dict, nucl_dict, aa_dict)

#%%
# write the clades TSV to be parsed by augur clades
write_clades_tsv(final_mutations_list)

#%%
# run augur clades to determine y_predictions 
run_augur_clades()


#%%
# compute the balanced accuracy of each clade's performance 
compute_performance(meta)


#%%
for i in final_mutations:
	print(f"{i}   :   {len(final_mutations[i])}  :  {len(clade_dict[i])}" )

#%%
for i in final_mutations:
	if i.startswith('2.3.4'):
		print(f"{i}   :   {len(final_mutations[i])}  :  {len(clade_dict[i])}" )





#%%

#%%
jsondict = json.load(open('../nextclade/clades/clades.json'))
predictions = [(id,jsondict['nodes'][id]['clade_membership']) for id in jsondict['nodes']]
predictions = pd.DataFrame(predictions, columns=['name','pred'])
predictions = predictions.loc[~predictions['name'].str.startswith("NODE")]

##%%
predictions = predictions.merge(meta, on='name')
N = predictions.shape[0]
## %% calculate metrics
check = (predictions['pred'] == predictions['clade'])
##%% perfect match 
print(check.sum() / N)
##%% calling a parent clade properly, but not the perfect match 
partial = predictions[['pred', 'clade']].apply(lambda x: x['clade'].startswith(x['pred']), axis=1)
print(partial.sum() / N)
## %%

wrong = predictions.loc[~check]
wrong = wrong.sort_values(['pred','clade'])
print(wrong[['pred','clade']].value_counts())
#%%
print(wrong['clade'].value_counts())
print(wrong['pred'].value_counts())




#%% ERROR ANALYSIS 
importlib.reload(augur_utils)
clade_designations = augur_utils.read_in_clade_definitions('../nextclade/clades/h5_clades.tsv')  # subtracts 1

tree = Phylo.read(input_tree, 'newick')
node_data = augur_utils.read_node_data([
	aa_json_path, 
	ancestral_json_path], input_tree)
all_muts = node_data['nodes']
mainref = augur_utils.get_reference_sequence_from_root_node(all_muts, tree.root.name)
# clade_membership = augur_utils.assign_clades(clade_designations, all_muts, tree, clade_dict['2.3.4.4'], mainref)

#%%
newtree = augur_utils.get_tree(clade_designations, all_muts, tree, mainref)

def fix_clade(tree, node_name):
	node = [x for x in tree.find_clades(order='preorder') if x.name == node_name][0]
	mutations = [[gene, site, alt] for gene in node.sequences for site,alt in node.sequences[gene].items()]
	for i in range(len(mutations)):
		mutations[i][1] += 1
	return [tuple(x) for x in mutations]
#%%
#final_mutations['2.3.4.4e'] = set(final_mutations['2.3.4.4e']).intersection(fix_clade(newtree, 'NODE_0005577'))
#final_mutations['2.2.2.1'] = fix_clade(newtree, 'NODE_0002423')
final_mutations['2.1.3.2'] = fix_clade(newtree, 'NODE_0001980')
#final_mutations['2.1.3.2a'] = fix_clade(newtree, 'NODE_0002186')
final_mutations['2.1.3.2a'] = fix_clade(newtree, 'NODE_0002284')
final_mutations['2.1.3.2b'] = fix_clade(newtree, 'NODE_0002335')

final_mutations['2.3.4.4'] = list(set(final_mutations['2.3.4.4']).intersection(fix_clade(newtree, 'NODE_0008237')))
# final_mutations['2.3.4.4e'] = fix_clade(newtree, 'NODE_0009126')
final_mutations['2.3.4.4e'] = fix_clade(newtree, 'NODE_0009214')
final_mutations['2.3.4.4d'] = fix_clade(newtree, 'NODE_0008521')
# final_mutations['2.3.4.4g'] = fix_clade(newtree, 'NODE_0008085')
final_mutations['2.3.4.4g'] = fix_clade(newtree, 'NODE_0008178')
final_mutations['1'] = list(set(final_mutations['1']).intersection(fix_clade(newtree, 'NODE_0001236')))
final_mutations['2.2.2'] = fix_clade(newtree, 'NODE_0002900')
final_mutations['2.2.2.1'] = fix_clade(newtree, 'NODE_0002964')