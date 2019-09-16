#!/bin/bash
set -e

fasta=$1
nreads=$2
ncpus=$3

full_fasta=$(realpath $fasta)
prefix=$(basename $fasta .fasta)

mkdir $prefix

# simulate reads
cd $prefix
# mason_simulator --num-threads $ncpus --seq-technology illumina -ir $full_fasta -n $nreads -o "${prefix}_left.fq" -or "${prefix}_right.fq" --seed $RANDOM
art_illumina -ss HS25 -na -i $full_fasta -p -l 150 -f 20 -m 200 -s 10 -o $prefix

# assemble
spades.py --threads $ncpus -1 "${prefix}1.fq" -2 "${prefix}2.fq" -o ./

# clean up reads
rm *.fq
rm -r K21
rm -r K33
rm -r K55
rm -r misc
rm -r corrected
rm -r tmp
rm assembly_graph.fastg
rm assembly_graph_with_scaffolds.gfa
rm before_rr.fasta
rm dataset.info
rm input_dataset.yaml
rm scaffolds.paths

# run prokka
# prokka_prefix="prokka_${prefix}"
# prokka --noanno --cpus $ncpus --outdir $prokka_prefix --prefix $prokka_prefix scaffolds.fasta

cd ..
