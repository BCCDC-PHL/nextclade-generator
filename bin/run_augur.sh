seqs=$1
meta=$2
threads=$3
refseq_name="NC_007362"

if [ $# -ne 3 ] && [ $# -ne 4 ]; then
	echo "USAGE: ./run_augur.sh [SEQUENCES_FASTA] [METADATA_TSV] [NTHREADS] (OUTDIR)" && 
	exit 0
fi 

if [ $# -eq 4 ]; then
	echo "Running in specified output directory." 
	base=$4
else
	echo "Running in internal nextclade directory." 
	base='../nextclade'
fi

mkdir -p $base/ancestral
mkdir -p $base/auspice
mkdir -p $base/clades
mkdir -p $base/inputs
mkdir -p $base/refine
mkdir -p $base/translate

echo "Starting sequence alignment ..." &&
augur align --sequences $seqs --nthreads $threads --reference-name ${refseq_name} --output $base/inputs/h5_align.fasta &&
printf "Completed sequence alignment. \nStarting tree construction ..." &&

augur tree --alignment $base/inputs/h5_align.fasta --nthreads $threads --output $base/inputs/h5_tree.nwk
printf "Completed tree construction. \nStarting augur refine ..." &&

augur refine --alignment $base/inputs/h5_align.fasta --tree $base/inputs/h5_tree.nwk --metadata $meta --output-tree $base/refine/h5_refined.nwk --output-node-data $base/refine/h5_node_data.json --root ${refseq_name} 

printf "Completed Augur refine. \nStarting Augur ancestral ..." &&

augur ancestral --tree $base/refine/h5_refined.nwk --alignment $base/inputs/h5_align.fasta --output-node-data $base/ancestral/h5_anc_node_data.json --output-sequences $base/ancestral/h5_ancestral.fasta 

printf "Completed Augur ancestral. \nStarting Augur translate ..." &&

augur translate --tree $base/refine/h5_refined.nwk --ancestral-sequences $base/ancestral/h5_anc_node_data.json --output-node-data $base/translate/h5_aa_node_data.json --alignment-output $base/translate/h5_align_%GENE.fasta --reference-sequence ../assets/h5n1_ref.gff3 

printf "Completed Augur translate." 

# augur clades --tree $base/refine/h5_refined.nwk --mutations $base/translate/h5_aa_node_data.json --reference ../assets/h5_ref.fasta --output-node-data $base/clades/h5_clades.json
