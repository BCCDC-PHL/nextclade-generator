import sys
from Bio import Phylo
import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
from itertools import islice
sys.path.append('/home/john.palmer/miniconda/envs/nextstrain/lib/python3.10/site-packages')
from augur.utils import read_node_data

def get_reference_sequence_from_root_node(all_muts, root_name):
	# attach sequences to root
	ref = {}
	try:
		ref['nuc'] = list(all_muts[root_name]["sequence"])
	except:
		print("WARNING in augur.clades: nucleotide mutation json does not contain full sequences for the root node.")

	if "aa_muts" in all_muts[root_name]:
		try:
			ref.update({gene:list(seq) for gene, seq in all_muts[root_name]["aa_sequences"].items()})
		except:
			print("WARNING in augur.clades: amino acid mutation json does not contain full sequences for the root node.")

	return ref

def read_in_clade_definitions(clade_file):
    '''
    Reads in tab-seperated file that defines clades by amino acid or nucleotide mutations
    Inheritance is allowed, but needs to be acyclic. Alleles can be overwritten by inheriting clades.
    Sites are 1 indexed in the file, and are converted to 0 indexed in the output
    Empty lines are ignored, comments after # are ignored
    Format::
        clade      gene    site     alt
        Clade_1    ctpE    81       D
        Clade_2    nuc     30642    T
        Clade_3    nuc     444296   A
        Clade_3    S       1        P
        # Clade_4 inherits from Clade_3
        Clade_4    clade   Clade_3
        Clade_4    pks8    634      T
        # Inherited allele can be overwritten
        Clade_4    S       1        L
    Parameters
    ----------
    clade_file : str
        meta data file
    Returns
    -------
    dict
        clade definitions as :code:`{clade_name:[(gene, site, allele),...]}`
    '''

    clades = defaultdict(lambda: defaultdict(str))
    df = pd.read_csv(
        clade_file,
        sep='\t' if clade_file.endswith('.tsv') else ',',
        comment='#'
    )

    clade_inheritance_rows = df[df['gene'] == 'clade']

    # Identify clades that inherit more than once
    clades_with_multiple_inheritance = clade_inheritance_rows[clade_inheritance_rows.duplicated(subset=["clade"])]['clade'].tolist()
    if len(clades_with_multiple_inheritance) > 0:
        raise ValueError(f"Clades {clades_with_multiple_inheritance} have multiple inheritance, that's not allowed")

    # Identify clades that inherit from non-existent clades
    missing_parent_clades = set(clade_inheritance_rows['site']) - set(df["clade"])
    if len(missing_parent_clades) > 0:
        raise ValueError(f"Clades {missing_parent_clades} are inherited from but are not defined")


    G = nx.DiGraph()

    # Use integer 0 as root so as not to conflict with any string clade names
    # String '0' can still be used this way
    root = 0
    # For every clade, add edge from root as default
    # This way all clades can be reached by traversal
    for clade in df.clade.unique():
        G.add_edge(root, clade)

    # Build inheritance graph
    # For clades that inherit, disconnect from root
    # Add edge from parent
    for _, row in clade_inheritance_rows.iterrows():
        G.remove_edge(root, row.clade)
        G.add_edge(row.site, row.clade)

    if not nx.is_directed_acyclic_graph(G):
        raise ValueError(f"Clade definitions contain cycles {list(nx.simple_cycles(G))}")

    # Traverse graph top down, so that children can inherit from parents and grandparents
    # Topological sort ensures parents are visited before children
    # islice is used to skip the root node (which has no parent)
    for clade in islice(nx.topological_sort(G),1,None):
        # Get name of parent clade
        # G.predecessors(clade) returns iterator, thus next() necessary
        # despite the fact that there should only be one parent
        parent_clade = next(G.predecessors(clade))
        # Inheritance from parents happens here
        # Allele dict is initialized with alleles from parent
        clades[clade] = clades[parent_clade].copy()
        for _, row in df[(df.clade == clade) & (df.gene != 'clade')].iterrows():
            # Overwrite of parent alleles is possible and happens here
            clades[clade][(row.gene, int(row.site)-1)] = row.alt

    # Convert items from dict[str, dict[(str,int),str]] to dict[str, list[(str,int,str)]]
    clades = {
        clade: [
            gene_site + (alt,)
            for gene_site, alt in clade_definition.items()
        ]
        for clade, clade_definition in clades.items()
        # If clause avoids root (helper) from being emmitted
        if clade != root
    }

    return clades


def is_node_in_clade(clade_alleles, node, ref):
	'''
	Determines whether a node matches the clade definition based on sequence
	For any condition, will first look in mutations stored in node.sequences,
	then check whether a reference sequence is available, and other reports 'non-match'
	Parameters
	----------
	clade_alleles : list
		list of clade defining alleles
	node : Bio.Phylo.BaseTree.Clade
		node to check, assuming sequences (as mutations) are attached to node
	ref : str or list
		positions
	Returns
	-------
	bool
		True if in clade
	'''
	conditions = []
	for gene, pos, clade_state in clade_alleles:
		pos = int(pos)
		if gene in node.sequences and pos in node.sequences[gene]:
			state = node.sequences[gene][pos]
		elif ref and gene in ref:
			state = ref[gene][pos]
		else:
			state = ''

		conditions.append(state==clade_state)

	return all(conditions)

def assign_clades(clade_designations, all_muts, tree, target_set, ref=None):
	'''
	Ensures all nodes have an entry (or auspice doesn't display nicely), tests each node
	to see if it's the first member of a clade (assigns 'clade_annotation'), and sets
	all nodes's clade_membership to the value of their parent. This will change if later found to be
	the first member of a clade.
	Parameters
	----------
	clade_designations :     dict
		clade definitions as :code:`{clade_name:[(gene, site, allele),...]}`
	all_muts : dict
		mutations in each node
	tree : Bio.Phylo.BaseTree.Tree
		phylogenetic tree to process
	ref : str or list, optional
		reference sequence to look up state when not mutated
	Returns
	-------
	dict
		mapping of node to clades
	'''

	clade_membership = {}

	# first pass to set all nodes to unassigned as precaution to ensure attribute is set
	for node in tree.find_clades(order = 'preorder'):
		clade_membership[node.name] = {'clade_membership': 'unassigned'}

	# count leaves
	for node in tree.find_clades(order = 'postorder'):
		node.leaf_count = 1 if node.is_terminal() else np.sum([c.leaf_count for c in node])

	for node in tree.get_nonterminals():
		for c in node:
			c.up=node
	tree.root.up = None
	tree.root.sequences = {'nuc':{}}
	tree.root.sequences.update({gene:{} for gene in all_muts[tree.root.name]['aa_muts']})

	# attach sequences to all nodes
	for node in tree.find_clades(order='preorder'):
		if node.up:
			node.sequences = {gene:muts.copy() for gene, muts in node.up.sequences.items()}
		for mut in all_muts[node.name]['muts']:
			a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]
			node.sequences['nuc'][pos] = d
		if 'aa_muts' in all_muts[node.name]:
			for gene in all_muts[node.name]['aa_muts']:
				for mut in all_muts[node.name]['aa_muts'][gene]:
					a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]

					if gene not in node.sequences:
						node.sequences[gene]={}
					node.sequences[gene][pos] = d

	#print(clade_membership['NODE_0005280'])
	# second pass to assign 'clade_annotation' to basal nodes within each clade
	# if multiple nodes match, assign annotation to largest
	# otherwise occasional unwanted cousin nodes get assigned the annotation
	print("Second pass.")
	for clade_name, clade_alleles in clade_designations.items():
		node_counts = []
		for node in tree.find_clades(order = 'preorder'):
			if node.name == 'NODE_0005577' and clade_name == '2.3.4.4e':
				print(node.name)
				print(sorted(node.sequences['nuc'].items()))
				print(sorted([x for x in clade_alleles if x[0] == 'nuc']))
			if is_node_in_clade(clade_alleles, node, ref):
				node_counts.append(node)
		# searching for the node with the largest number of leaves
		sorted_nodes = sorted(node_counts, key=lambda x: x.leaf_count, reverse=True)
		if clade_name == '2.3.4.4e':
			print(sorted_nodes)
		if len(sorted_nodes) > 0:
			target_node = sorted_nodes[0]

			clade_membership[target_node.name] = {'clade_annotation': clade_name, 'clade_membership': clade_name}
	
	print(clade_membership['NODE_0005577'])
	# third pass to propagate 'clade_membership'
	# don't propagate if encountering 'clade_annotation'
	for node in tree.find_clades(order = 'preorder'):
		for child in node:
			if 'clade_annotation' not in clade_membership[child.name]:
				#print(child.name)
				# if child.name == 'EPI_ISL_295501':
				# 	print(f'Child = EPI_ISL_295501')
				# 	print(f'Parent = {node.name}')
					
				clade_membership[child.name]['clade_membership'] = clade_membership[node.name]['clade_membership']
	print(clade_membership['NODE_0005577'])
	return clade_membership

def get_tree(clade_designations, all_muts, tree, ref=None):

	clade_membership = {}

	# first pass to set all nodes to unassigned as precaution to ensure attribute is set
	for node in tree.find_clades(order = 'preorder'):
		clade_membership[node.name] = {'clade_membership': 'unassigned'}

	# count leaves
	for node in tree.find_clades(order = 'postorder'):
		node.leaf_count = 1 if node.is_terminal() else np.sum([c.leaf_count for c in node])

	for node in tree.get_nonterminals():
		for c in node:
			c.up=node
	tree.root.up = None
	tree.root.sequences = {'nuc':{}}
	tree.root.sequences.update({gene:{} for gene in all_muts[tree.root.name]['aa_muts']})

	# attach sequences to all nodes
	for node in tree.find_clades(order='preorder'):
		if node.up:
			node.sequences = {gene:muts.copy() for gene, muts in node.up.sequences.items()}
		for mut in all_muts[node.name]['muts']:
			a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]
			node.sequences['nuc'][pos] = d
		if 'aa_muts' in all_muts[node.name]:
			for gene in all_muts[node.name]['aa_muts']:
				for mut in all_muts[node.name]['aa_muts'][gene]:
					a, pos, d = mut[0], int(mut[1:-1])-1, mut[-1]

					if gene not in node.sequences:
						node.sequences[gene]={}
					node.sequences[gene][pos] = d

	# second pass to assign 'clade_annotation' to basal nodes within each clade
	# if multiple nodes match, assign annotation to largest
	# otherwise occasional unwanted cousin nodes get assigned the annotation
	print("Second pass.")
	for clade_name, clade_alleles in clade_designations.items():
		node_counts = []
		for node in tree.find_clades(order = 'preorder'):
			if is_node_in_clade(clade_alleles, node, ref):
				node_counts.append(node)
		# searching for the node with the largest number of leaves
		sorted_nodes = sorted(node_counts, key=lambda x: x.leaf_count, reverse=True)
		if len(sorted_nodes) > 0:
			target_node = sorted_nodes[0]

			clade_membership[target_node.name] = {'clade_annotation': clade_name, 'clade_membership': clade_name}
	
	# third pass to propagate 'clade_membership'
	# don't propagate if encountering 'clade_annotation'
	for node in tree.find_clades(order = 'preorder'):
		for child in node:
			if 'clade_annotation' not in clade_membership[child.name]:
				#print(child.name)
				# if child.name == 'EPI_ISL_295501':
				# 	print(f'Child = EPI_ISL_295501')
				# 	print(f'Parent = {node.name}')
					
				clade_membership[child.name]['clade_membership'] = clade_membership[node.name]['clade_membership']
	return tree