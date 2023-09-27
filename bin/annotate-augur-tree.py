#%%
import json 
import os, sys
import pandas as pd 
import numpy as np
import math
from collections import Counter
import re 

EXPR_HIERARCHICAL = re.compile('^[\d\.]+[a-z]?$')

# %% helper functions 
def print_leaves(node):
	if 'children' not in node or not node['children']:
		print(node)
		return 
	
	for child in node.get('children', []):
		print_leaves(child)


	

	
def longest_common_prefix(strlist):
	prefix = ''
	strlist = sorted(strlist)
	first = strlist[0]
	last = strlist[-1]

	for i in range(min(len(first), len(last))):
		if first[i] != last[i]:
			return prefix
		prefix += first[i]
	
	return prefix

#%%
def assign_clade(node, clade):
	if 'node_attrs' not in node:
		node['node_attrs'] = {}
	if 'clade_membership' not in node['node_attrs']:
		node['node_attrs']['clade_membership'] = {}
	node['node_attrs']['clade_membership']['value'] = clade

def assign_clade_leaves(node, clade_dict):
	if 'children' not in node or not node['children']:
		assign_clade(node, clade_dict.get(node['name'],'unassigned'))
		return 
	
	for child in node.get('children', []):
		assign_clade_leaves(child, clade_dict)

#%%
def assign_clade_nodes(node):
	if 'children' not in node or not node['children']:
		clade = node['node_attrs']['clade_membership'].get('value','')
		return [clade]
	

	clades = []
	for child in node.get('children', []):
		c = assign_clade_nodes(child)
		clades += c

	

	# next, try the cumulative summation method 
	# threshold = int(N * 0.75)
	# clade_counter = Counter(clades)

	# cumsum = 0
	# clade_subset = []
	# # finds the most essential clades such that their cumulative sum totals over 75% 
	# for clade, count in clade_counter.most_common():
	# 	cumsum += count
	# 	clade_subset += [clade]
	# 	if cumsum > threshold:
	# 		break


	N = len(clades)
	threshold = 0.05


	clade_counter = Counter(clades)
	print(clade_counter)
	clade_subset = [x for x, y in clade_counter.items() if y / N > threshold]
	

	# try to compute again on the subset of most prevalent clades
	prefix = longest_common_prefix(clade_subset)

	if node['name'] == 'NODE_0012729':
		print(node['name'], prefix, clade_subset)

	if prefix != '':
		assign_clade(node, prefix)
		return clades

	# prefix, _ = clade_counter.most_common(1)[0]
	# if node['name'] == 'NODE_0012729':
	# 	print(node['name'], prefix, clade_counter)

	# assign the clade 
	# assign_clade(node, prefix)
	
	return clades 

#%%
def count_subtree(node, node_name, found):
	if 'children' not in node or not node['children']:
		return 

	if node['name'] == node_name:
		f = True
	elif found == True:
		f = True
	else:
		f = False

	if f:
		print(node['name'])

	for child in node.get('children', []):
		count_subtree(child, node_name, f)

	return 


#%%
def get_top_node_per_clade(node, node_dict):
	if 'children' not in node or not node['children']:
		return 

	for child in node.get('children', []):
		get_top_node_per_clade(child, node_dict)

	if 'clade_membership' in node['node_attrs']:
		clade = node['node_attrs']['clade_membership']['value']
		node_dict[clade] = node['name']
	return 


def extract_clade_mutations(node, target_nodes, clade_mutations):
	if 'children' not in node or not node['children']:
		return 
	
	if node['name'] in target_nodes:
		clade_mutations[target_nodes[node['name']]] = node['branch_attrs']['mutations']

	for child in node.get('children', []):
		extract_clade_mutations(child, target_nodes, clade_mutations)

	return 

#%% main operations 
base_path = '/home/john.palmer/git/nextclade-generator/nextclade/'

treepath = '/home/john.palmer/work/flu/ref/augur_h5n1/ha_dataset/master2/auspice/tree.json'
treepath = base_path + 'auspice/tree.json'
with open(treepath, 'r') as handle:
	data = json.load(handle)
tree = data['tree']

#%%
metapath = base_path + 'inputs/master-meta-final-nodup.tsv'
df = pd.read_csv(metapath, sep='\t').set_index('name')
df['clade'] = df['clade'].replace({np.nan:''})
clade_dict = df['clade'].to_dict()
clade_counts = df['clade'].value_counts().to_dict()

#%%
assign_clade_leaves(tree, clade_dict)


#%% output tree into new JSON file 
# data['tree'] = tree

# with open(base_path + 'auspice/tree_tip_label.json', 'w') as outfile:
# 	json.dump(data, outfile)

#%%
assign_clade_nodes(tree)

#%%
clade_top_nodes = {}
get_top_node_per_clade(tree, clade_top_nodes)

clade_top_nodes = {y:x for x, y in clade_top_nodes.items()}
#%%
clade_mutations = {}
extract_clade_mutations(tree, clade_top_nodes, clade_mutations)

#%%

# %% first pass to assign clades to all tips 

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade

def json_to_clade(json_node):
	"""Recursively convert a JSON tree node to a Biopython Clade."""
	name = json_node.get('name')
	branch_length = json_node.get('branch_length', None)
	clade = Clade(name=name, branch_length=branch_length)

	children = json_node.get('children', [])
	for child_json in children:
		child_clade = json_to_clade(child_json)
		clade.clades.append(child_clade)

	return clade

def json_to_phylo(json_tree):
	"""Convert a JSON tree to a Biopython Phylo tree."""
	clade = json_to_clade(json_tree)
	tree = Phylo.BaseTree.Tree.from_clade(clade)
	return tree

# Example usage:
tree_json = {
	"name": "root",
	"branch_length": 0,
	"children": [
		{
			"name": "A",
			"branch_length": 1,
			"children": [
				{"name": "A1", "branch_length": 1},
				{"name": "A2", "branch_length": 1}
			]
		},
		{
			"name": "B",
			"branch_length": 1,
			"children": [
				{"name": "B1", "branch_length": 1},
				{"name": "B2", "branch_length": 1}
			]
		}
	]
}

phylo_tree = json_to_phylo(tree_json)

# Display the tree
Phylo.draw_ascii(phylo_tree)
# %%
