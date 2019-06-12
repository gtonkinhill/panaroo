# Output

### gene_presence_absence.csv

A csv file describing which gene is in which sample. If a gene cluster is present in a sample, the sequence name of the representative for that sample is given in the matrix. The corresponding DNA and protein sequence can then be matched to those found in the `combined_DNA_CDS.fasta` and `combined_protein_CDS.fasta` files. The format is the same as that given by [Roary](https://sanger-pathogens.github.io/Roary/).

### gene_mobility.csv

A csv file indicating the degree (number of edges) of each gene in the graph. This indicates the number of different genetic contexts this gene has been identified in. A degree of two suggests that the gene is only seen in one context.

### struct_presence_absence.csv

A csv file which lists the presence and abscence of different genomic rearrangement events. The genes involved in each event are listed in the respective column names of the csv. The thresholds for calling these events can be changed by adjusting the `--min_edge_support_sv` parameter when calling panaurus.

### pan_genome_reference.fa

This is a similar output to that produced by [Roary](https://sanger-pathogens.github.io/Roary/). It creates a linear reference genome of all the genes found in the dataset. The order of the genes in this reference are not significant.

### final_graph.gml

The final pan-genome graph generated by Panaurus. This can be viewed easily using Cytoscape (see XX for more details). It includes all the metainformation such as gene annotation and which gene/edge is present in which genome. It can also be loaded into python using networkx for further processing.

### gene_data.csv

This is a very large file mainly used internally in the program. It links each gene sequnece and annotation to the internal representations used. It can be useful in interpreting some of the output especially the 'final_graph.gml' file.

### combined_DNA_CDS.fasta

This is a fasta file which includes all nucleotide sequence for both the annotated genes and those refound by the program. The gene names are the internal ones used by Panaurus. These can be translated to the original names using the 'gene_data.csv' file.

### combined_protein_CDS.fasta

Similar to the `combined_DNA_CDS.fasta` file, this is a fasta file which includes all protein sequence for both the annotated genes and those refound by the program. The gene names are the internal ones used by Panaurus. These can be translated to the original names using the 'gene_data.csv' file.