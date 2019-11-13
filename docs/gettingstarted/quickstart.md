# Quick Start


## Basic usage

Using GFFs in the same format as output by Prokka run:

```
mkdir results
panaroo -i *.gff -o results
```

## Mode

By default Panaroo runs in its strictest (most conservative) mode. We have found that for most use cases this removes potential sources of contamination and error whilst retaining the majority of genes researchers are interested in. 

Very rare plasmids are difficult to distguish from contamination. Thus, if you are interested in retaining such plasmids at the expense of added contamination we recommend running panaroo using its most sensitive mode

```
panaroo -i *.gff -o results --mode relaxed
```


## Pre-processing scripts

It is usually a good idea to perform some rudimentary quality checks on the input data prior to running Panaroo. This can be accomplished using the pre-processing script provided.

The reference mash database can be downloaded from https://mash.readthedocs.io/en/latest/tutorials.html

```
panaroo-qc -t 3 --graph_type all -i *.gff --ref_db refseq.genomes.k21.s1000.msh -o results
```


## Running Prokka

Whilst Panaroo is designed to correct for many sources of error introduced during annotation, you can improve the consistency of annotations by using the same gene model in Prokka. We have provided a helpful script to do this:

```
mkdir reannotated
run_prokka -i *.fasta -o reannotated
```

