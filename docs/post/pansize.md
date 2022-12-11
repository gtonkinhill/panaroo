# Estimating the Pangenome Size

Panaroo includes script for estimating the pangenome size as well as gene gain and loss rates according to both the Finitely Many Genes (FMG) model and the Infinitely Many Genes (IMG) model. These models are described in a number of papers listed in the references at the bottom of this page.

These model based approaches are preferable to the common accumulation curves often used in pangenome analyses. Such curves fail to account for the phylogenetic relationship between isolates and thus it is inappropriate to use them to compare different data sets.

### Infinitely Many Genes model

The IMG model allows for gene gain from an unbounded reservoir of new genes, and gene loss. As the reservoir is unbounded, the same gene gain event can only occur once. This might represent horizontal gene transfer from a diverged taxa with gene loss representing the conversion of genes to pseudogenes or deletion in reproduction. This model is described in Baumdicker et al. 2012 and Collins et al. 2012.

To estimate the parameters of this model, a dated phylogeny based on the core genome is required. Such phylogenies can be produced using [BEAST](https://www.beast2.org/) or by combining faster methods such as [IQ-TREE](http://www.iqtree.org/) and [BactDating](https://xavierdidelot.github.io/BactDating/)

An implementation following that in Collins et al. 2012 is given in Panaroo and can be run as

```
mkdir img_results
panaroo-img --pa gene_presence_absence.Rtab -o img_results -tree dated_phylogeny.newick
```

#### Parameters

For a full description of the parameters and model please see the Collins et al. 2012 paper.

```
Estimate model parameters for either the Infinitely Many Genes Model using
gene frequencies

optional arguments:
  -h, --help            show this help message and exit
  --tree TREE           A dated phylogeny.
  --pa PRESENCE_ABSENCE
                        A presence/absence produced by Panaroo.
  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
                        location of an output directory
  -D {1,2}              Number of separate rate classes to use for the
                        dispensable genome. Can be either 1 or 2.
  --no_essential        Removes essential gene class from model
  --no_constraint       Removes constraint that u/v must equal the genome
                        size.
  --model {coalescent,fixed}
                        Model to fit. Can be one of 'coalescent' or 'fixed'.
  --fit {cp,gf}         Whether to use the gene frequency spectrum or the
                        core/pangenome size for fitting (default=gf)
  --init_u1 U1          Initial value for u1 (default=0.01).
  --init_u2 U2          Initial value for u2 (default=0.01).
  --init_v1 V1          Initial value for v1 (default=0.01).
  --init_v2 V2          Initial value for v2 (default=0.01).
  --init_ess GESS       Initial value for the number of essential genes
                        (default=2000).
  --verbose             print additional output
  --version             show program's version number and exit
```


### Finitely Many Genes model

The FMG model differs from the IMG mode in that the same gene can be gained and lost multiple times. Moreover, the pool of genes from which a new gene can be drawn is finite. The FMG model is a version of the standard two-state phylogenetic models which are often used to study the evolution of binary traits.

We have implemented a simplified version of that described in Zamani-Dahaj et al. 2016.

Similar to the IMG model a dated phylogeny based on the core genome is required. The script can then be run as

```
panaroo-fmg --tree dated_phylogeny.newick --pa gene_presence_absence_renamed.Rtab -o fmg_results.txt
```

#### Parameters

For a full description of the paramters and model please see the Zamani-Dahaj et al. 2016 paper.

```
Estimate model parameters for the Finitely Many Genes Model

optional arguments:
  -h, --help            show this help message and exit
  --tree TREE           A dated phylogeny.
  --pa PRESENCE_ABSENCE
                        A presence/absence produced by Panaroo.
  -o OUTPUTFILE, --outfile OUTPUTFILE
                        Name of outputfile.
  --nboot NBOOT         The number of sub-sampling bootstrap iterations to
                        perform.
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --verbose             print additional output
  --version             show program's version number and exit
```

### Alternative Tools

The output from Panaroo can also be used as input to a number of other post-processing tools for estimating gene gain and loss rates. These include

- [PanX](http://pangenome.de/)
- [panicmage](http://www.baumdickerlab.de/index.php/software/panicmage)
- [pantagruel](https://github.com/flass/pantagruel)


### References

Collins,R.E. and Higgs,P.G. (2012) Testing the Infinitely Many Genes Model for the Evolution of the Bacterial Core Genome and Pangenome. Mol. Biol. Evol., 29, 3413–3425.

Baumdicker,F., Hess,W.R. and Pfaffelhuber,P. (2012) The infinitely many genes model for the distributed genome of bacteria. Genome Biol. Evol., 4, 443–456.

Zamani-Dahaj,S.A., Okasha,M., Kosakowski,J. and Higgs,P.G. (2016) Estimating the Frequency of Horizontal Gene Transfer Using Phylogenetic Models of Gene Gain and Loss. Mol. Biol. Evol., 33, 1843–1857.
