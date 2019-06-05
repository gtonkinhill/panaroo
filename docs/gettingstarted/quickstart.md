# Quick Start


## Basic usage

Using GFFs in the same format as output by Prokka run:

```
mkdir results
panaurus -i *.gff -o results
```

### Running Prokka

You can make use of our handy script to run Prokka using the same prodigal model on each genome. This helps to give more consistent results between genomes:

```
mkdir reannotated
run_prokka -i *.fasta -o reannotated
```

