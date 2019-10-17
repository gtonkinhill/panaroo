#!/bin/bash
set -e

#conda activate panaroo

fasta=$1
nreads=$2
ncpus=$3

full_fasta=$(readlink -f $fasta)
prefix=$(basename $fasta .fasta)

cd ngs_art_sim

mkdir $prefix

# simulate reads
cd $prefix
# mason_simulator --num-threads $ncpus --seq-technology illumina -ir $full_fasta -n $nreads -o "${prefix}_left.fq" -or "${prefix}_right.fq" --seed $RANDOM
# art_illumina -ss HS25 -na -i $full_fasta -p -l 150 -f 20 -m 200 -s 10 -o $prefix
art_illumina -ss HS25 -i $full_fasta -o ${prefix}_ -l 150 -f 20 -p -m 200 -s 10 --noALN

# assemble
spades.py --phred-offset 33 --threads $ncpus -1 "${prefix}_1.fq" -2 "${prefix}_2.fq" -o ./

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

cd ../..


# run prokka
prokka --proteins pan_sim_*prokka_DB.fasta --cpus $ncpus --outdir ./prokka_art_assem/${prefix} --prefix $prefix ./ngs_art_sim/${prefix}/scaffolds.fasta
