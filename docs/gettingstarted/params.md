# Parameters

#### Mode

The main parameter that should be adjusted when using Panaroo is the `mode`. This changes a number of defaults to either favour `strict`, `moderate` or `sensitive` filtering.
The default parameter settings for each of these modes is given in the table at the end of this page. 

#### Gene Alignment

Panaroo can generate multiple sequence alignments from the resulting gene clusters. Currently, the user can choose from either [PRANK](https://www.ncbi.nlm.nih.gov/pubmed/24170401),  [Clustal](http://www.clustal.org/omega/) or the default, [MAFFT](https://mafft.cbrc.jp/alignment/software/). MSA can be run either for the core genome using `-a core` or for every gene cluster using `-a pan`. The frequency of a gene in the your sample required to classify it as 'core' can be set using `--core_threshold`.

Thus to align all genes present in at least 98% of isolates using clustal and 10 cpus you would run Panaroo as

```
panaroo -i *.gff -o ./results/ --clean-mode strict -a core --aligner clustal --core_threshold 0.98 -t 10
```

You can also output unaligned gene sequences by specifying `--aligner none`. Additionally, user @revinci has provided a separate script for generating alignments after running Panaroo, which is described [here](https://github.com/gtonkinhill/panaroo/issues/306).

#### Cluster Thresholds

The Panaroo algorithm initially performs a conservative clustering step before collapsing genes into possible families. It is usually best to use the dafault parameters for this initial clustering stage.

Thus we recommend using the defaults for `--threshold` (0.98) and `--len_dif_percent` (0.98).

If you wish to adjust the level at which Panaroo collapses genes into putative families we suggest changing the family sequence identity level (default 0.7). Thus to run Panaroo using a more relaxed threshold of 50% identity you could run

```
panaroo -i *.gff -o ./results/ --clean-mode strict -f 0.5
```

#### Paralogs

Panaroo splits paralogs into separate clusters by default. Merging paralogs can be enabled by running Panaroo as

```
panaroo -i *.gff -o ./results/  --clean-mode strict --merge_paralogs
```

#### Refinding Genes

In order to identify genes that have been missed by annotation software, Panaroo incorporates a refinding step. Suppose two clusters geneA and geneB are adjacent in the Panaroo pangenome graph. If geneA is present in a genome but its neighbour (geneB) is not then Panaroo searches the sequence surrounding geneA for the presence of geneB. The radius of this search in nucleotides is controlled by `--search_radius`, with the default being 5000.  

As such missing genes are often the results of assembly fragmentation, the refinding step only requires that a proportion of the missing gene is located. This proportion can be controlled using the `--refind_prop_match` parameter.

Thus to refind genes that match at least 50% of the sequence within a radius of 1000 nucleotides you could run

```
panaroo -i *.gff -o ./results/ --clean-mode strict --refind_prop_match 0.5 --search_radius 1000
```


### Panaroo Usage

```
usage: panaroo [-h] -i INPUT_FILES [INPUT_FILES ...] -o OUTPUT_DIR
               --clean-mode {strict,moderate,sensitive}
               [--remove-invalid-genes] [-c ID] [-f FAMILY_THRESHOLD]
               [--len_dif_percent LEN_DIF_PERCENT]
               [--family_len_dif_percent FAMILY_LEN_DIF_PERCENT]
               [--merge_paralogs] [--search_radius SEARCH_RADIUS]
               [--refind_prop_match REFIND_PROP_MATCH]
               [--refind-mode {default,strict,off}]
               [--min_trailing_support MIN_TRAILING_SUPPORT]
               [--trailing_recursive TRAILING_RECURSIVE]
               [--edge_support_threshold EDGE_SUPPORT_THRESHOLD]
               [--length_outlier_support_proportion LENGTH_OUTLIER_SUPPORT_PROPORTION]
               [--remove_by_consensus {True,False}]
               [--high_var_flag CYCLE_THRESHOLD_MIN]
               [--min_edge_support_sv MIN_EDGE_SUPPORT_SV]
               [--all_seq_in_graph] [--no_clean_edges] [-a {core,pan}]
               [--aligner {prank,clustal,mafft,none}] [--codons]
               [--core_threshold CORE] [--core_subset SUBSET]
               [--core_entropy_filter HC_THRESHOLD] [-t N_CPU]
               [--codon-table TABLE] [--quiet] [--version]

panaroo: an updated pipeline for pangenome investigation

optional arguments:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --codon-table TABLE   the codon table to use for translation (default=11)
  --quiet               suppress additional output
  --version             show program's version number and exit

Input/output:
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        input GFF3 files (usually output from running Prokka).
                        Can also take a file listing each gff file line by
                        line.
  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
                        location of an output directory

Mode:
  --clean-mode {strict,moderate,sensitive}
                        The stringency mode at which to run panaroo. Must be
                        one of 'strict','moderate' or 'sensitive'. Each of
                        these modes can be fine tuned using the additional
                        parameters in the 'Graph correction' section.

                        strict:
                        Requires fairly strong evidence (present in  at least
                        5% of genomes) to keep likely contaminant genes. Will
                        remove genes that are refound more often than they were
                        called originally.

                        moderate:
                        Requires moderate evidence (present in  at least 1% of
                        genomes) to keep likely contaminant genes. Keeps genes
                        that are refound more often than they were called
                        originally.

                        sensitive:
                        Does not delete any genes and only performes merge and
                        refinding operations. Useful if rare plasmids are of
                        interest as these are often hard to disguish from
                        contamination. Results will likely include  higher
                        number of spurious annotations.
  --remove-invalid-genes
                        removes annotations that do not conform to the
                        expected Prokka format such as those including
                        premature stop codons.

Matching:
  -c ID, --threshold ID
                        sequence identity threshold (default=0.98)
  -f FAMILY_THRESHOLD, --family_threshold FAMILY_THRESHOLD
                        protein family sequence identity threshold
                        (default=0.7)
  --len_dif_percent LEN_DIF_PERCENT
                        length difference cutoff (default=0.98)
  --family_len_dif_percent FAMILY_LEN_DIF_PERCENT
                        length difference cutoff at the gene family level
                        (default=0.0)
  --merge_paralogs      don't split paralogs

Refind:
  --search_radius SEARCH_RADIUS
                        the distance in nucleotides surronding the neighbour
                        of an accessory gene in which to search for it
  --refind_prop_match REFIND_PROP_MATCH
                        the proportion of an accessory gene that must be found
                        in order to consider it a match
  --refind-mode {default,strict,off}
                        The stringency mode at which to re-find genes.

                        default:
                        Will re-find similar gene sequences. Allows for
                        premature stop codons and incorrect lengths to account
                        for misassemblies.

                        strict:
                        Prevents fragmented, misassembled, or potential
                        pseudogene sequences from being re-found.

                        off:
                        Turns off all re-finding steps.

Graph correction:
  --min_trailing_support MIN_TRAILING_SUPPORT
                        minimum cluster size to keep a gene called at the end
                        of a contig
  --trailing_recursive TRAILING_RECURSIVE
                        number of times to perform recursive trimming of low
                        support nodes near the end of contigs
  --edge_support_threshold EDGE_SUPPORT_THRESHOLD
                        minimum support required to keep an edge that has been
                        flagged as a possible mis-assembly
  --length_outlier_support_proportion LENGTH_OUTLIER_SUPPORT_PROPORTION
                        proportion of genomes supporting a gene with a length
                        more than 1.5x outside the interquatile range for
                        genes in the same cluster (default=0.01). Genes
                        failing this test will be re-annotated at the shorter
                        length
  --remove_by_consensus {True,False}
                        if a gene is called in the same region with similar
                        sequence a minority of the time, remove it. One of
                        'True' or 'False'
  --high_var_flag CYCLE_THRESHOLD_MIN
                        minimum number of nested cycles to call a highly
                        variable gene region (default = 5).
  --min_edge_support_sv MIN_EDGE_SUPPORT_SV
                        minimum edge support required to call structural
                        variants in the presence/absence sv file
  --all_seq_in_graph    Retains all DNA sequence for each gene cluster in the
                        graph output. Off by default as it uses a large amount
                        of space.
  --no_clean_edges      Turn off edge filtering in the final output graph.

Gene alignment:
  -a {core,pan}, --alignment {core,pan}
                        Output alignments of core genes or all genes. Options
                        are 'core' and 'pan'. Default: 'None'
  --aligner {prank,clustal,mafft,none}
                        Specify an aligner. Options:'prank', 'clustal', and
                        default: 'mafft'
  --codons              Generate codon alignments by aligning sequences at the
                        protein level
  --core_threshold CORE
                        Core-genome sample threshold (default=0.95)
  --core_subset SUBSET  Randomly subset the core genome to these many genes
                        (default=all)
  --core_entropy_filter HC_THRESHOLD
                        Manually set the Block Mapping and Gathering with
                        Entropy (BMGE) filter. Can be between 0.0 and 1.0. By
                        default this is set using the Tukey outlier method.
```

#### Default Parameters

$N$ is the number of genomes in the dataset.


| Parameter                | Description                                                                                           | strict                   | moderate                 | relaxed |
|--------------------------|-------------------------------------------------------------------------------------------------------|--------------------------|--------------------------|---------|
| min_trailing_support   | Minimum number of genomes supporting a trailing cluster.  Cluster with less than this will be deleted | $\max(2,\lceil 0.05 N \rceil)$ | $\max(2,\lceil 0.01 N \rceil)$ | 2       |
| trailing_recursive     | Maximum number of times to recursively delete trailing clusters with degree $1$ that have less than the minimum trailing support. | $\infty$                   | $\infty$                   | 1       |
| min_edge_support_sv    | Minimum edge support required to be included in `struct_presence_absence.csv` | $\max(2,\lceil 0.01 N\rceil)$ | $\max(2,\lceil 0.01 N\rceil)$ | 2       |
| remove_by_consensus    | Whether or not to remove clusters that are refound more times than they were originally annotated. This can indicate a poor initial annotation in a minority of genomes. | True                     | False                    | False   |
| edge_support_threshold | This only affects the GML output. Edges that are below this threshold and which join clusters of higher support are removed. This improves the readability of the graph. | $\max(2,\lceil 0.01 N\rceil)$ | $\max(2,\lceil 0.01 N\rceil)$ | 0       |
