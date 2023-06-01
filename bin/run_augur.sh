seqs=$1
threads=$2

mkdir -p out_align

augur align --sequences $seqs --nthreads $threads --reference-name NC_007362 --output out_align/h5_align.fasta

mkdir -p out_tree
augur tree --alignment out_align/h5_align.fasta --nthreads $threads --output out_tree/h5_tree.nwk

