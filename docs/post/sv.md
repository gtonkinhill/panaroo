# Structural Variant Analysis

As Panaroo builds a full representation of each gene in the pangenome we can use this to investigate gene rearrangements. 

To simplify association analyses between these rearrangments and phenotypes of interest we concentrate on gene triplets. As paralogs are split in the pangenome graph each sample can only be present in one path. Thus each genome can only contain one version of triplet with the same central gene. This is illustrated in the figure below.

We record the presence and absence of these gene triplets in each sample in the `struct_presence_absence.csv` file. This can be used as input to [Pyseer](https://pyseer.readthedocs.io/en/master/) in a similar way to the gene presence/absence matrix.

<img src="_figures/structural_variants.png" width="800">

### Example

We continue on from the example given in the [Gene Association Study](/post/pwas.md) section.

```
for anti in AZM     CRO     CFM     CIP     PEN     SMX     TET
do
python ~/pyseer/pyseer-runner.py --lmm --phenotypes ./metadata/pathogenwatch-neigo-eurogasp2013-amr-profile-1.tab --pres ./panaroo_output/struct_presence_absence.Rtab --similarity ./pyseer_out/phylogeny_K.tsv --phenotype-column $anti --output-patterns ./pyseer_out/struct_patterns_${anti}.txt > ./pyseer_out/${anti}_struct_was.txt
done
```

again counting the number of patterns for false discovery analysis

```
python ~/pyseer/scripts/count_patterns.py ./pyseer_out/struct_patterns_AZM.txt
```

which gives

```
Patterns:       3239
Threshold:      1.54E-05
```

We can now summarise the rearrangement results using *R*

```{r}
gono_struct_was <- do.call(rbind, map(antibiotics, function(ant){
  tbl <- fread(paste(c("./gono_harris_gwas/",ant,"_struct_was.txt"), collapse = ""), data.table = FALSE)
  tbl$antibiotic <- ant
  return(tbl)
}))

gono_struct_was <- gono_struct_was[order(gono_struct_was$`lrt-pvalue`),]
gono_struct_was <- gono_struct_was[!grepl("bad-chisq", gono_struct_was$notes),]

# threshold form running count_patterns in pyseer
sig_threshold <- 0.05/(3239*length(antibiotics))

sum(gono_struct_was$`lrt-pvalue`<sig_threshold)
top_results <- gono_struct_was[gono_struct_was$`lrt-pvalue`<sig_threshold,]
top_results <- top_results[!duplicated(top_results$variant),]

write.csv(top_results, file="./gono_harris_gwas/gono_struct_was_annotated.csv", quote = FALSE)
```
