mkdir -p refine
mkdir -p ancestral
mkdir -p translate
mkdir -p clades

meta=$1

augur refine --alignment out_align/h5_align.fasta --tree out_tree/h5_tree.nwk --metadata $meta --output-tree refine/h5_refined.nwk --output-node-data refine/h5_node_data.json --root NC_007362 

augur ancestral --tree refine/h5_refined.nwk --alignment out_align/h5_align.fasta --output-node-data ancestral/h5_anc_node_data.json --output-sequences ancestral/h5_ancestral.fasta 

augur translate --tree refine/h5_refined.nwk --ancestral-sequences ancestral/h5_anc_node_data.json --output-node-data translate/h5_aa_node_data.json --alignment-output translate/h5_align_%GENE.fasta --reference-sequence h5_ref.gff 

# augur clades --tree ../refine/h5_refined.nwk --mutations ../translate/h5_aa_node_data.json --reference reference.fasta --output-node-data h5_clades.json
