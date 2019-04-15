
###################################################################
Command line arguments
###################################################################


Panaurus
========

-h, --help            show this help message and exit
-c ID, --threshold ID
                      sequence identity threshold (default=0.95)
-f FAMILY_THRESHOLD, --family_threshold FAMILY_THRESHOLD
                      protein family sequence identity threshold
                      (default=0.7)
--len_dif_percent LEN_DIF_PERCENT
                      length difference cutoff (default=0.95)
-i INPUT_FILES, --input INPUT_FILES
                      input GFF3 files (usually output from running Prokka)
-o OUTPUT_DIR, --out_dir OUTPUT_DIR
                      location of an output directory
--min_trailing_support MIN_TRAILING_SUPPORT
                      minimum cluster size to keep a gene called at the end
                      of a contig (default=2)
--trailing_recursive TRAILING_RECURSIVE
                      number of times to perform recursive triming of low
                      support nodes near the end of contigs (default=2)
--max_cycle_size MAX_CYCLE_SIZE
                      maximum cycle size for collapsing gene families
                      (default=20)
--min_edge_support_sv MIN_EDGE_SUPPORT_SV
                      minimum edge support required to call structural
                      variants in the presence/absence sv file
--no_split            don't split paralogs
--mode
                      the stringency mode at which to run panaurus. One of
                      'strict', 'moderate' or 'relaxed' (default='strict')
-t N_CPU, --threads N_CPU
                      number of threads to use (default=1)
--verbose             print additional output



run_prokka
=========

-h, --help            show this help message and exit
-i INPUT_FILES, --input INPUT_FILES
                      input GFF3 files (usually output from running Prokka)
-o OUTPUT_DIR, --out_dir OUTPUT_DIR
                      location of an output directory
--add_prokka_cmds ADD_PROKKA_CMDS
                      additional commands to supply to Prokka (these are not
                      checked!)
--num_training NUM_TRAINING
                      number of genomes to use in training prodigal
                      (default=10)
-t N_CPU, --threads N_CPU
                      number of threads to use (default=1)
--force               overwrite old commands
--verbose             print additional output
