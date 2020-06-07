#!/usr/bin/env python3
"""
A script to merge a single GFF file with a pre-existing panaroo output
"""

import os
import tempfile
import networkx as nx
import shutil
import sys

from .__init__ import __version__
from .prokka import process_prokka_input
from .cdhit import run_cdhit
from .generate_network import generate_network
from .isvalid import *
from .merge_graphs import *

def get_options(args): #options for integrating (combination of merge graph and cdhit options
   
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
        type=str)
    
    io_opts.add_argument(
        "-i",
        "--input_gff",
        dest="input_gff",
        required=True,
        help="input gff file of new genome to be integrated",
        type=str)
        
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of a new output directory",
                         type=lambda x: is_valid_folder(parser, x))

    matching = parser.add_argument_group('Matching')

    matching.add_argument("-c",
                          "--threshold",
                          dest="id",
                          help="sequence identity threshold (default=0.95)",
                          default=0.98,
                          type=float)

    matching.add_argument(
        "-f",
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        default=0.7,
        type=float)

    matching.add_argument("--len_dif_percent",
                          dest="len_dif_percent",
                          help="length difference cutoff (default=0.95)",
                          default=0.95,
                          type=float)

    matching.add_argument("--merge_paralogs",
                          dest="merge_paralogs",
                          help="don't split paralogs",
                          action='store_true',
                          default=False)

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
        default=0.1)

    parser.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=(
            "minimum edge support required to call structural variants" +
            " in the presence/absence csv file (default=max(2, 0.01*n_samples))"
        ),
        default=2,
        type=int)

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
        default=None)
    
    core.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'prank', 'clustal', and default: 'mafft'",
        type=str,
        choices=['prank', 'clustal', 'mafft'],
        default="mafft")
    
    core.add_argument("--core_threshold",
                      dest="core",
                      help="Core-genome sample threshold (default=0.95)",
                      type=float,
                      default=0.95)
    
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
                        default=1)
    
    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False)
                        
    parser.add_argument("--dirty",
                          dest="dirty",
                          help="keep temporary directory containing cluster files and cdhit output",
                          action='store_true',
                          default=False)
                          
    parser.add_argument("--version'",
                        action="version",
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
   
    return (args)

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

def reformat_network(single_gml, output_dir, isolateName): #Generate network output needs to be reformatted for single gff inputs
    
    """Reformats the output of generate_network() for linear graphs to allow input into merge_graphs()"""
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
        
        to_replace = {"[": "", "]": "", "'": ""}
        y['centroid'] = replace_all(str(y['centroid']),to_replace)
        y['dna'] = replace_all(str(y['dna']),to_replace)
        y['protein'] = replace_all(str(y['protein']),to_replace)
        
        y['hasEnd'] = int(y['hasEnd'])
        y['mergedDNA'] = int(y['mergedDNA'])
        y['paralog'] = int(y['paralog'])
        
        y['longCentroidID'] = list(y['longCentroidID'])
        y['seqIDs'] = list(y['seqIDs'])
    
    single_gml.graph.update({'isolateNames' : 'x'}) # isolateName from gff filename
    nx.write_gml(single_gml, output_dir + "final_graph.gml")
    
    return single_gml

def merge_graphs(directories, filename, n_cpu, temp_dir, len_dif_percent, Id, family_threshold, length_outlier_support_proportion, quiet, merge_paralogs, output_dir, min_edge_support_sv, aln, alr, core):

   """Merge graphs in directories"""
    print(
        "Merging graphs is still under active development and may change frequently!"
    )
    # Load graphs
    print("Loading graphs...")
    graphs, isolate_names, id_mapping = load_graphs(
        [d + "final_graph.gml" for d in directories], n_cpu=n_cpu)

    # cluster centroids
    print("Clustering centroids...")
    clusters, seqid_to_centroid, orig_ids, ids_len_stop = cluster_centroids(
        graphs=graphs,
        outdir=temp_dir,
        directories=directories,
        id_mapping=id_mapping,
        len_dif_percent=len_dif_percent,
        identity_threshold=Id,
        n_cpu=n_cpu)

    # perform initial merge
    print("Performing inital merge...")
    G = simple_merge_graphs(graphs, clusters)

    print("Number of nodes in merged graph: ", G.number_of_nodes())

    print("Collapsing at DNA...")
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          dna_error_threshold=0.98,
                          correct_mistranslations=True,
                          length_outlier_support_proportion=
                          length_outlier_support_proportion,
                          n_cpu=n_cpu,
                          quiet=quiet)[0]

    print("Number of nodes in merged graph: ", G.number_of_nodes())

    print("Collapsing at families...")
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          family_threshold=family_threshold,
                          correct_mistranslations=False,
                          length_outlier_support_proportion= length_outlier_support_proportion,
                          n_cpu=n_cpu,
                          quiet=quiet)[0]

    print("Number of nodes in merged graph: ", G.number_of_nodes())

    # Generate output
    print("Generating output...")

    # if requested merge paralogs
    if merge_paralogs:
        G = merge_paralogs(G)

    G.graph['isolateNames'] = isolate_names
    mems_to_isolates = {}
    for i, iso in enumerate(isolate_names):
        mems_to_isolates[i] = iso

    if not quiet:
        print("writing output...")

    # write out roary like gene_presence_absence.csv
    # get original annotaiton IDs, lengths and whether or
    # not an internal stop codon is present
    orig_ids = {}
    ids_len_stop = {}
    for i, d in enumerate(directories):
        with open(d + "gene_data.csv", 'r') as infile:
            next(infile)
            for line in infile:
                line = line.split(",")
                if line[2] not in id_mapping[i]: continue  #its been filtered
                orig_ids[id_mapping[i][line[2]]] = line[3]
                ids_len_stop[id_mapping[i][line[2]]] = (len(line[4]),
                                                        "*" in line[4][1:-3])

    G = generate_roary_gene_presence_absence(G,
                                             mems_to_isolates=mems_to_isolates,
                                             orig_ids=orig_ids,
                                             ids_len_stop=ids_len_stop,
                                             output_dir=output_dir)
    #Write out presence_absence summary
    generate_summary_stats(output_dir=output_dir)

    # write pan genome reference fasta file
    generate_pan_genome_reference(G,
                                  output_dir=output_dir,
                                  split_paralogs=False)

    # write out common structural differences in a matrix format
    generate_common_struct_presence_absence(
        G,
        output_dir=output_dir,
        mems_to_isolates=mems_to_isolates,
        min_variant_support=min_edge_support_sv)

    # add helpful attributes and write out graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['centroid'] = ";".join(G.nodes[node]['centroid'])
        G.nodes[node]['dna'] = ";".join(conv_list(G.nodes[node]['dna']))
        G.nodes[node]['protein'] = ";".join(conv_list(
            G.nodes[node]['protein']))
        G.nodes[node]['genomeIDs'] = ";".join(
            [str(m) for m in G.nodes[node]['members']])
        G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
        G.nodes[node]['degrees'] = G.degree[node]
        G.nodes[node]['members'] = list(G.nodes[node]['members'])
        G.nodes[node]['seqIDs'] = list(G.nodes[node]['seqIDs'])

    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            [str(m) for m in G.edges[edge[0], edge[1]]['members']])
        G.edges[edge[0],
                edge[1]]['members'] = list(G.edges[edge[0],
                                                   edge[1]]['members'])
    
    for index, name in enumerate(G.graph['isolateNames']): #Corrects isolate name for single gff being returned as list
        if name == 'x':
            G.graph['isolateNames'][index] = filename
            
    nx.write_gml(G, output_dir + "final_graph.gml")

    # write out merged gene_data and combined_DNA_CDS files
    with open(output_dir + "gene_data.csv", 'w') as outdata, \
    open(output_dir + "combined_DNA_CDS.fasta", 'w') as outdna:
        for i, d in enumerate(directories):
            with open(d + "gene_data.csv", 'r') as infile:
                header = next(infile)
                if i == 0: outdata.write(header)
                for line in infile:
                    line = line.strip().split(",")
                    if line[2] not in id_mapping[i]:
                        continue  #its been filtered
                    line[2] = id_mapping[i][line[2]]
                    outdata.write(",".join(line) + "\n")
                    outdna.write(">" + line[2] + "\n")
                    outdna.write(line[5] + "\n")

    # #Write out core/pan-genome alignments
    if aln == "pan":
        if not quiet: print("generating pan genome MSAs...")
        generate_pan_genome_alignment(G, temp_dir, output_dir, n_cpu,
                                      alr, isolate_names)
        core_nodes = get_core_gene_nodes(G, core, len(isolate_names))
        concatenate_core_genome_alignments(core_nodes, output_dir)
    elif aln == "core":
        print('yes')
        if not quiet: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, output_dir,
                                       n_cpu, alr, isolate_names,
                                       core, len(isolate_names))
    return 

def main(): #Takes a single GFF input, generates a graph and merges with a pre-existing graph
    args = get_options(sys.argv[1:])
    
    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
      os.mkdir(args.output_dir)
         
    args.input_dir = os.path.join(args.input_dir, "")
    args.output_dir = os.path.join(args.output_dir, "")
        
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    
    directories = [args.input_dir, temp_dir]
    
    gff_file = [args.input_gff]
 
    filename = os.path.basename(args.input_gff).split(".")[0]
    
    print("Processing input")
    process_prokka_input(gff_list = gff_file,
                         output_dir = temp_dir,
                         quiet = args.quiet, 
                         n_cpu = args.n_cpu) 
    
    cd_hit_out = temp_dir + "combined_protein_cdhit_out.txt"
    
    run_cdhit(input_file=temp_dir + "combined_protein_CDS.fasta",
              output_file=cd_hit_out,
              id=args.id,
              quiet=args.quiet,
              n_cpu=args.n_cpu)
    
    print("Generating network")
    
    single_gml, centroid_contexts_single, seqid_to_centroid_single = generate_network(cluster_file=cd_hit_out + ".clstr",
                                                                                    data_file=temp_dir + "gene_data.csv",
                                                                                    prot_seq_file=temp_dir + "combined_protein_CDS.fasta", 
                                                                                    all_dna=args.all_seq_in_graph)
    print("Reformatting network")
    reformat_network(single_gml = single_gml,
                     output_dir = temp_dir,
                     isolateName = filename)
    
    merge_graphs(directories, 
                 filename, 
                 args.n_cpu, 
                 temp_dir, 
                 args.len_dif_percent,
                 args.id,
                 args.family_threshold, 
                 args.length_outlier_support_proportion, 
                 args.quiet,
                 args.merge_paralogs, 
                 args.output_dir, 
                 args.min_edge_support_sv, 
                 args.aln,
                 args.alr, 
                 args.core)

    #remove temporary directory if dirty = True
    if not args.dirty:
        shutil.rmtree(temp_dir)
    
    sys.exit(0)

if __name__ == '__main__':
    main()
