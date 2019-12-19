# Pangenome Association Study

Panaroo produces output that can be provided to [Pyseer](https://pyseer.readthedocs.io/en/master/) to perform pangenome association analyses. To demonstrate this we provide an example analysis of antibiotic resistance in *Neisseria gonorrhoeae*.

We look at the association between gene presence/absence and antibiotic resistance in a large [European sample](https://doi.org/10.1016/S1473-3099(18)30225-1) of 1054 isolates. Data required to replicate the analysis is downloadable at [Pathogenwatch (eurogasp2013)](https://pathogen.watch/collection/eurogasp2013).


### Annotation

The assemblies from the study were first annotated using Prokka

```
mkdir prokka_output
for fasta in ./genomes/*.fasta
do
prefix=$(basename ${fasta} .fasta)
prokka --cpus 5 --genus Neisseria --usegenus --outdir ./prokka_output/${prefix} --prefix $prefix $fasta
done
```

### Panaroo

We can now run Panaroo, setting the option to build a multiple sequence alignment of the core genes using MAFFT.

```
mkdir panaroo_output
panaroo -i ./prokka_output/*/*.gff -o panaroo_output -t 24 --verbose -a core
```

In order to control for population structure in our association analyses we build a phylogeny from the core gene alignment using [Iqtree](http://www.iqtree.org/). 

```
cd panaroo_output
iqtree -s core_gene_alignment.aln -pre core_tree -nt 24 -fast -m GTR
cd ..
```

### Pyseer

We are now ready to run [Pyseer](https://pyseer.readthedocs.io/en/master/). We first look at association between antibiotic resistance and gene presence/absence.

```
python ~/pyseer/scripts/phylogeny_distance.py --lmm ./panaroo_output/core_tree.treefile > pyseer_out/phylogeny_K.tsv

for anti in AZM     CRO     CFM     CIP     PEN     SMX     TET
do
python ~/pyseer/pyseer-runner.py --lmm --phenotypes ./metadata/pathogenwatch-neigo-eurogasp2013-amr-profile-1.tab --pres ./panaroo_output/gene_presence_absence.Rtab --similarity ./pyseer_out/phylogeny_K.tsv --phenotype-column $anti --output-patterns ./pyseer_out/gene_patterns_${anti}.txt > ./pyseer_out/${anti}_gwas.txt
done
```

We also have to count the number of patterns in order to control for multiple testing.

```
python ~/pyseer/scripts/count_patterns.py ./pyseer_out/gene_patterns_AZM.txt
```

This results in:

```
Patterns:       529
Threshold:      9.45E-05
```

We now have everything we need to start looking at some results! We make use of a small bit of *R* code here to summarise the results.

```{r}
antibiotics <- c('AZM', 'CRO', 'CFM', 'CIP', 'PEN', 'SMX', 'TET')

gono_gwas <- do.call(rbind, map(antibiotics, function(ant){
  tbl <- fread(paste(c("./gono_harris_gwas/",ant,"_gwas.txt"), collapse = ""), data.table = FALSE)
  tbl$antibiotic <- ant
  return(tbl)
}))

gono_gwas <- gono_gwas[order(gono_gwas$`lrt-pvalue`),]
gono_gwas <- gono_gwas[!grepl("bad-chisq", gono_gwas$notes),]

# threshold form running count_patterns in pyseer
sig_threshold <- 0.05/(529*length(antibiotics))

sum(gono_gwas$`lrt-pvalue`<sig_threshold)
sig_hits <- gono_gwas[gono_gwas$`lrt-pvalue`<sig_threshold,]

write.csv(sig_hits, file = "pan_was_results.csv", quote = FALSE, row.names = FALSE)
```

