# Quality Control

Panaroo provides a pre-processing script which helps to identify possible issues with samples prior to running the full pangenome algorithm. Such checks can also be useful for many other analysis tasks.

The script can be called using 3 cpus as 

```
mkdir results
panaroo-qc -i *.gff -o results -t 3 --graph_type all --ref_db refseq.genomes.k21.s1000.msh
```

### Parameters

```
Generate quality control plots prior to a Panaroo run

optional arguments:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --graph_type {contam,ngenes,all,mds,ncontigs}
                        the type of graph to generate (default='all')
  --ref_db REF_DB       reference mash database for contamination
                        quantification.
  --version             show program's version number and exit

Input/output:
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        input GFF3 files (usually output from running Prokka)
  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
                        location of an output directory
```

The following section includes examples from running the algorithm on a large [Klebsiella](https://www.pnas.org/content/112/27/E3574) dataset. 

### Contamination

The script wraps the popular [Mash](https://mash.readthedocs.io/en/latest/tutorials.html) algorithm, and produces plots indicating possible sources of contamination as well as outliers.


### Multi-Dimensional Scaling (MDS)

After using Mash to create a pairwise distance matrix between each sample, the script generates an MDS plot. Samples that appear as outliers in this plot should be investigated for possible contamination, poor assembly or other problems.

<iframe seamless="seamless" scrolling="yes" src="_figures/MDS_mash_plot.html" width="800" height="500" frameborder="0"></iframe>

### Outliers

The script also produces plots which help to highlight samples with an unusual number of contigs or genes. This can indicate a problem with the sample.

<iframe seamless="seamless" scrolling="yes" src="_figures/ngenes_boxplot.html" width="800" height="500" frameborder="0"></iframe>


