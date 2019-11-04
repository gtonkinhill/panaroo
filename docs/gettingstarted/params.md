# Parameters

```
usage: panaroo [-h] -i INPUT_FILES [INPUT_FILES ...] -o OUTPUT_DIR [-c ID]
               [-f FAMILY_THRESHOLD] [--len_dif_percent LEN_DIF_PERCENT]
               [--merge_paralogs] [--search_radius SEARCH_RADIUS]
               [--refind_prop_match REFIND_PROP_MATCH]
               [--mode {strict,moderate,relaxed}]
               [--min_trailing_support MIN_TRAILING_SUPPORT]
               [--trailing_recursive TRAILING_RECURSIVE]
               [--edge_support_threshold EDGE_SUPPORT_THRESHOLD]
               [--remove_by_consensus {True,False}]
               [--high_var_flag CYCLE_THRESHOLD_MIN]
               [--min_edge_support_sv MIN_EDGE_SUPPORT_SV]
               [--all_seq_in_graph] [-a ALN] [--aligner ALR]
               [--core_threshold CORE] [-t N_CPU] [--verbose] [--version]

panaroo: an updated pipeline for pan-genome investigation

optional arguments:
  -h, --help            show this help message and exit
  -t N_CPU, --threads N_CPU
                        number of threads to use (default=1)
  --verbose             print additional output
  --version             show program's version number and exit

Input/output:
  -i INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
                        input GFF3 files (usually output from running Prokka)
  -o OUTPUT_DIR, --out_dir OUTPUT_DIR
                        location of an output directory

Matching:
  -c ID, --threshold ID
                        sequence identity threshold (default=0.95)
  -f FAMILY_THRESHOLD, --family_threshold FAMILY_THRESHOLD
                        protein family sequence identity threshold
                        (default=0.7)
  --len_dif_percent LEN_DIF_PERCENT
                        length difference cutoff (default=0.95)
  --merge_paralogs      don't split paralogs

Refind:
  --search_radius SEARCH_RADIUS
                        the distance in nucleotides surronding the neighbour
                        of an accessory gene in which to search for it
  --refind_prop_match REFIND_PROP_MATCH
                        the proportion of an accessory gene that must be found
                        in order to consider it a match

Graph correction:
  --mode {strict,moderate,relaxed}
                        the stringency mode at which to run panaroo. One of
                        'strict', 'moderate' or 'relaxed' (default='strict')
  --min_trailing_support MIN_TRAILING_SUPPORT
                        minimum cluster size to keep a gene called at the end
                        of a contig
  --trailing_recursive TRAILING_RECURSIVE
                        number of times to perform recursive triming of low
                        support nodes near the end of contigs
  --edge_support_threshold EDGE_SUPPORT_THRESHOLD
                        minimum support required to keep and edge that has
                        been flagged as a possible mis-assembly
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

Gene alignment:
  -a ALN, --alignment ALN
                        Output alignments of core genes or all genes. Options
                        are 'core' and 'pan'. Default: 'None'
  --aligner ALR         Specify an aligner. Options:'prank', 'clustal', and
                        default: 'mafft'
  --core_threshold CORE
                        Core-genome sample threshold (default=0.95)
```
