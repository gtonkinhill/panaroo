# Parameters

```
usage: panaroo [-h] -i INPUT_FILES [INPUT_FILES ...] -o OUTPUT_DIR [-c ID]
                [-f FAMILY_THRESHOLD] [--len_dif_percent LEN_DIF_PERCENT]
                [--merge_paralogs] [--mode {strict,moderate,relaxed}]
                [--min_trailing_support MIN_TRAILING_SUPPORT]
                [--trailing_recursive TRAILING_RECURSIVE]
                [--edge_support_diff EDGE_SUPPORT_DIFF]
                [--remove_by_consensus]
                [--min_edge_support_sv MIN_EDGE_SUPPORT_SV] [-a ALN]
                [--aligner ALR] [--core_threshold CORE] [-t N_CPU] [--verbose]
                [--version]

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
  --edge_support_diff EDGE_SUPPORT_DIFF
                        maximum fraction difference between an edge's support
                        and those of the nodes it connects
  --remove_by_consensus
                        if a gene is called in the same region with similar
                        sequence a minority of the time, remove it
  --min_edge_support_sv MIN_EDGE_SUPPORT_SV
                        minimum edge support required to call structural
                        variants in the presence/absence sv file

Gene alignment:
  -a ALN, --alignment ALN
                        Output alignments of core genes or all genes. Options
                        are 'core' and 'pan'. Default: 'None'
  --aligner ALR         Specify an aligner. Options:'prank', 'clustal', and
                        default: 'mafft'
  --core_threshold CORE
                        Core-genome sample threshold (default=0.95)
```
