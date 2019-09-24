# Quick Start


## Basic usage

Using GFFs in the same format as output by Prokka run:

```
mkdir results
panaroo -i *.gff -o results
```

## Pre-processing scripts

It is usually a good idea to perform some rudimentary quality checks on the input data prior to running Panaroo. This can be accomplished using the pre-processing script provided.

The reference mash database can be downloaded from https://mash.readthedocs.io/en/latest/tutorials.html

```
python pan-qc-runner.py -t 3 --graph_type all -i *.gff --ref_db refseq.genomes.k21.s1000.msh -o results
```


### Running Prokka

You can make use of our handy script to run Prokka using the same prodigal model on each genome. This helps to give more consistent results between genomes:

```
mkdir reannotated
run_prokka -i *.fasta -o reannotated
```

