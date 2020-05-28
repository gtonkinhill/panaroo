#!/usr/bin/env python3
"""
A script to merge a single GFF file with a pre-existing panaroo output
"""

import os
import tempfile
import networkx as nx
import subprocess
import shutil
import numpy as np

from .__init__ import __version__
from .prokka import process_prokka_input
from .cdhit import run_cdhit
from .generate_network import generate_network
from .isvalid import *

def get_options(): #options for integrating (combination of merge graph and cdhit options
   
    import argparse

    description = 'Integrate new gff file into pre-existing graph'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo_integrate')

    io_opts = parser.add_argument_group('Input/output')
    
    io_opts.add_argument(
        "-d",
        "--input_dir",
        dest="input_dir",
        required=True,
        help="input directory for gml of pre-existing panaroo output",
        type=str) #SPECIFY THE LOCATIONS OF THE PREEXISTING GRAPH 
    
    io_opts.add_argument(
        "-i",
        "--input_gff",
        dest="input_gff",
        required=True,
        help="input gff file of new genome to be integrated",
        type=str) #SPECIFY THE GFF TO BE INTEGRATED  
        
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of a new output directory",
                         type=lambda x: is_valid_folder(parser, x)) #SPECIFY THE OUTPUT DIRECTORY OF THE FINAL GRAPH. USED IN PROCESS_PROKKA_INPUT, RUN_CDHIT, NX.WRITE_GML

    matching = parser.add_argument_group('Matching')

    matching.add_argument("-c",
                          "--threshold",
                          dest="id",
                          help="sequence identity threshold (default=0.95)",
                          default=0.98,
                          type=float) #SEQUENCE IDENTITY THRESHOLD. USED IN RUN-CDHIT

    matching.add_argument(
        "-f",
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        default=0.7,
        type=float) #PROTEIN SEQUENCE IDENTITY THRESHOLD. USED IN RUN-CDHIT-EST

    matching.add_argument("--len_dif_percent",
                          dest="len_dif_percent",
                          help="length difference cutoff (default=0.95)",
                          default=0.95,
                          type=float) #USED IN CLUSTER CENTROIDS AND RUN-CDHIT

    matching.add_argument("--merge_paralogs",
                          dest="merge_paralogs",
                          help="don't split paralogs",
                          action='store_true',
                          default=False) #USED IN MERGE_PARALOGS FROM CLEAN NETWORK  

    matching.add_argument(
        "--length_outlier_support_proportion",
        dest="length_outlier_support_proportion",
        help=
        ("proportion of genomes supporting a gene with a length more " +
         "than 1.5x outside the interquatile range for genes in the same cluster"
         +
         " (default=0.01). Genes failing this test will be re-annotated at the "
         + "shorter length"),
        type=float,
        default=0.1) #USED IN COLLAPSE FAMILIES IN CLEAN NETWORK  

    parser.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=(
            "minimum edge support required to call structural variants" +
            " in the presence/absence csv file (default=max(2, 0.01*n_samples))"
        ),
        default=2,
        type=int) #USED IN generate_common_struct_presence_absence IN GENERATE OUTPUT

    # MSA options
    core = parser.add_argument_group('Gene alignment')
    core.add_argument(
        "-a",
        "--alignment",
        dest="aln",
        help=("Output alignments of core genes or all genes. Options are" +
              " 'core' and 'pan'. Default: 'None'"),
        type=str,
        choices=['core', 'pan'],
        default=None) #USED IN generate_pan_genome_alignment AND generate_core_genome_alignment FROM GENERATE_OUTPUT
    
    core.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'prank', 'clustal', and default: 'mafft'",
        type=str,
        choices=['prank', 'clustal', 'mafft'],
        default="mafft") #USED IN generate_pan_genome_alignment AND generate_core_genome_alignment FROM GENERATE_OUTPUT
    
    core.add_argument("--core_threshold",
                      dest="core",
                      help="Core-genome sample threshold (default=0.95)",
                      type=float,
                      default=0.95) #USED IN generate_core_genome_alignment and get_core_gene_nodes FROM GENERATE_OUTPUT
    
    graph = parser.add_argument_group('Graph correction')
    
    graph.add_argument("--all_seq_in_graph",
                       dest="all_seq_in_graph",
                       help=("Retains all DNA sequence for each gene cluster in the graph " +
                             "output. Off by default as it uses a large amount of space."),
                       action='store_true',
                       default=False)
        
    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1) #USED IN LOAD_GRAPHS, cluster_centroids, collapse_families, generate_pan_genome_alignment, generate_core_genome_alignment,PROCESS_PROKKA INPUT, RUN_CDJIT
    
    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False) #USED IN collapse_families, Write out core/pan-genome alignments, PROCESS_PROKKA_INPUT 
        
    parser.add_argument("--version'",
                        action="version",
                        version='%(prog)s ' + __version__) #IMPORTED FROM INNIT.PY

    args = parser.parse_args()
    return (args)


def reformat_network(single_gml, output_dir, isolateName): #Generate network output needs to be reformatted for single gff inputs
    
    for adj in single_gml._adj:
        for x in single_gml._adj[adj]:
            y = single_gml._adj[adj][x]
        
            y.pop('members')
            zero = {'members': 0}
            y.update(zero)
        
            genomes = {'genomeIDs' : '0'}
            y.update(genomes)
        
    for node in single_gml._node:
        y = single_gml._node[node]
        y.pop('members')
        
        zero = {'members': 0} #members are assigned intbitset[0]. needs to be 0
        y.update(zero)
        
        rep = "[']"
        for char in rep:
            y['centroid'] = (str(y['centroid'])).replace(char, '')
            y['dna'] = (str(y['dna'])).replace(char, '')
            y['protein'] = (str(y['protein'])).replace(char, '')
        
        y['hasEnd'] = int(y['hasEnd'])
        y['mergedDNA'] = int(y['mergedDNA'])
        y['paralog'] = int(y['paralog'])
        
        y['longCentroidID'] = list(y['longCentroidID'])
        y['seqIDs'] = list(y['seqIDs'])
    
    single_gml.graph.update({'isolateNames' : 'x'}) # isolateName from gff filename
    nx.write_gml(single_gml, output_dir + "final_graph.gml")
    
    return single_gml


def main(): #Takes a single GFF input, generates a graph and merges with a pre-existing graph
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
     
    directories = args.input_dir + ' ' + temp_dir

    gff_file = []
    gff_file.append(args.input_gff)
            
    filename = (str(args.input_gff).split('/')[-1]).split('.')[0]
 
    process_prokka_input(gff_list = gff_file, 
                         output_dir = temp_dir, 
                         quiet = args.quiet, 
                         n_cpu = args.n_cpu) 
    
    cd_hit_out = temp_dir + "/combined_protein_cdhit_out.txt"
    
    run_cdhit(input_file=temp_dir + "/combined_protein_CDS.fasta",
              output_file=cd_hit_out,
              id=args.id,
              quiet=args.quiet,
              n_cpu=args.n_cpu)
    
    single_gml, centroid_contexts_single, seqid_to_centroid_single = generate_network(cluster_file=cd_hit_out + ".clstr",
                                                                                    data_file=temp_dir + "/gene_data.csv",
                                                                                    prot_seq_file=temp_dir + "/combined_protein_CDS.fasta",
                                                                                    all_dna=args.all_seq_in_graph)
    
    reformat_network(single_gml = single_gml,
                     output_dir = temp_dir,
                     isolateName = filename)
    
    merge_command = 'panaroo-merge -d ' + directories
    merge_command += ' -o ' + args.output_dir
    merge_command += ' -f ' + str(args.family_threshold)
    merge_command += ' --len_dif_percent ' + str(args.len_dif_percent)
    
    if args.merge_paralogs == True:
        merge_command += ' ' + str(args.merge_paralogs)
    merge_command += ' --length_outlier_support_proportion ' + str(args.length_outlier_support_proportion)
    merge_command += ' --min_edge_support_sv ' + str(args.min_edge_support_sv)
    if args.aln is not None:
        merge_command += ' -a ' + str(args.aln)
    merge_command += ' --aligner ' + args.alr
    merge_command += ' --core_threshold ' + str(args.core)
    merge_command += ' -t ' + str(args.n_cpu)
    if args.quiet == True:
        merge_command += ' ' + str(args.quiet)
    #merge_command += ' ' + str(args.version)
    
    subprocess.run(merge_command, shell = True)
    
    merged = nx.read_gml(args.output_dir + "/final_graph.gml")

    for index, name in enumerate(merged.graph['isolateNames']): #Corrects isolate name for single gff being returned as list
        if name == 'x':
            merged.graph['isolateNames'][index] = filename

    nx.write_gml(merged, args.output_dir + "/final_graph.gml")
    
    #remove temporary directory
    shutil.rmtree(temp_dir)
    
    return

if __name__ == '__main__':
    main()
