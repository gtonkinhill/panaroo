# Quick Start


## Basic usage

Using GFFs in the same format as output by Prokka run:

```
mkdir results
panaroo -i *.gff -o results --clean-mode strict
```

If you are using GFFs from RefSeq they can occasionally include annotations that do not conform the expected. This is usually due to a premature stop codon or a gene of invalid length. By default Panaroo will fail to parse these annotations. However, you can set Panaroo to ignore invalid annotaion by enabling the `remove-invalid-genes` flag 

```
panaroo -i *.gff -o results --clean-mode strict --remove-invalid-genes
```

## Mode

By default Panaroo runs in its strictest (most conservative) mode. We have found that for most use cases this removes potential sources of contamination and error whilst retaining the majority of genes researchers are interested in. 

Very rare plasmids are difficult to distguish from contamination. Thus, if you are interested in retaining such plasmids at the expense of added contamination we recommend running panaroo using its most sensitive mode

```
panaroo -i *.gff -o results --clean-mode sensitive
```


## Pre-processing scripts

It is usually a good idea to perform some rudimentary quality checks on the input data prior to running Panaroo. This can be accomplished using the pre-processing script provided.

The reference mash database can be downloaded from https://mash.readthedocs.io/en/latest/tutorials.html

```
panaroo-qc -t 3 --graph_type all -i *.gff --ref_db refseq.genomes.k21s1000.msh -o results
```


## Running Prokka

Whilst Panaroo is designed to correct for many sources of error introduced during annotation, you can improve the consistency of annotations by using the same gene model in Prokka. We have provided a helpful script to do this:

```
mkdir reannotated
run_prokka -i *.fasta -o reannotated
```

