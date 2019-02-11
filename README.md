# squeaky
An updated pipeline for pan-genome investigation


### Pipeline

1. Run prodigal on full dataset (so as not to run into issue with its model inference step). We will also take Prokka output as input.
2. Cluster using cd-hit with option -G 0
3. Build a pan-population pan-genome graph using the cd-hit clusters and adjaceny information from the assemblies
4. Split paralogs into multiple nodes
5. Trim genes that have low support and appear at the ends of the graph. These are likely false positive occuring as a result of the difficulty in calling genes near the end of contigs.
6. Collapse 'bubbles' into gene families based on a second more relaxed cut-off. Perhaps use fastbaps here.
7. Prepare suitable output files for further downstream analysis
    * Core genome alignment (split out as a seperate script)
    * Gene presence/absence (or count) matrix with and without genes being collapsed into families
    * Pan genome reference fasta
    * Optional MSA of each gene cluster (split out as a seperate script)
    * Output for visulaisation in Gephi (add paths if requested)

### Potential downstream analysis

1. Use capture-reacture methods that take into account measurement error (ghost records) http://www.math.canterbury.ac.nz/~r.vale/Mta.pdf or https://arxiv.org/pdf/1401.3290.pdf
2. Analyse the number of 'forks' in the graph which represent recombination events as we have already split out paralogs.
3. Possible identify phages
4. Include R code for the analysis of gene presence/absence.
5. Produce some function for common plots that people like.
