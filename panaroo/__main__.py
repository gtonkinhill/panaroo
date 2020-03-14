import os, sys
import tempfile
from Bio import SeqIO
import shutil
import networkx as nx
import ast

from .isvalid import *
from .set_default_args import set_default_args
from .prokka import process_prokka_input
from .cdhit import check_cdhit_version
from .cdhit import run_cdhit
from .generate_network import generate_network
from .generate_output import *
from .clean_network import *
from .find_missing import find_missing
from .generate_alignments import check_aligner_install

from .__init__ import __version__


def get_options(args):
    import argparse

    description = 'panaroo: an updated pipeline for pangenome investigation'
    parser = argparse.ArgumentParser(description=description, prog='panaroo')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help=("input GFF3 files (usually output from running Prokka). " +
              "Can also take a file listing each gff file line by line."),
        type=str,
        nargs='+')
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=lambda x: is_valid_folder(parser, x))

    matching = parser.add_argument_group('Matching')
    matching.add_argument("-c",
                          "--threshold",
                          dest="id",
                          help="sequence identity threshold (default=0.95)",
                          type=float)
    matching.add_argument(
        "-f",
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        type=float)
    matching.add_argument("--len_dif_percent",
                          dest="len_dif_percent",
                          help="length difference cutoff (default=0.95)",
                          type=float)
    matching.add_argument("--merge_paralogs",
                          dest="merge_paralogs",
                          help="don't split paralogs",
                          action='store_true',
                          default=False)

    refind = parser.add_argument_group('Refind')
    refind.add_argument(
        "--search_radius",
        dest="search_radius",
        help=("the distance in nucleotides surronding the " +
              "neighbour of an accessory gene in which to search for it"),
        default=5000,
        type=int)
    refind.add_argument(
        "--refind_prop_match",
        dest="refind_prop_match",
        help=("the proportion of an accessory gene that must " +
              "be found in order to consider it a match"),
        default=0.2,
        type=float)

    graph = parser.add_argument_group('Graph correction')
    graph.add_argument(
        "--mode",
        dest="mode",
        help=("the stringency mode at which to run panaroo. One of 'strict'" +
              ", 'moderate' or 'relaxed' (default='strict')"),
        choices=['strict', 'moderate', 'relaxed'],
        default='strict')
    graph.add_argument(
        "--min_trailing_support",
        dest="min_trailing_support",
        help=("minimum cluster size to keep a gene called at the " +
              "end of a contig"),
        type=int)
    graph.add_argument(
        "--trailing_recursive",
        dest="trailing_recursive",
        help=("number of times to perform recursive trimming of low support " +
              "nodes near the end of contigs"),
        type=int)
    graph.add_argument(
        "--edge_support_threshold",
        dest="edge_support_threshold",
        help=(
            "minimum support required to keep and edge that has been flagged" +
            " as a possible mis-assembly"),
        type=float)
    graph.add_argument(
        "--remove_by_consensus",
        dest="remove_by_consensus",
        type=ast.literal_eval,
        choices=[True, False],
        help=
        ("if a gene is called in the same region with similar sequence a minority "
         + "of the time, remove it. One of 'True' or 'False'"),
        default=None)
    graph.add_argument(
        "--high_var_flag",
        dest="cycle_threshold_min",
        help=(
            "minimum number of nested cycles to call a highly variable gene " +
            "region (default = 5)."),
        type=int,
        default=5)
    graph.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=("minimum edge support required to call structural variants" +
              " in the presence/absence sv file"),
        type=int)
    graph.add_argument(
        "--all_seq_in_graph",
        dest="all_seq_in_graph",
        help=("Retains all DNA sequence for each gene cluster in the graph " +
              "output. Off by default as it uses a large amount of space."),
        action='store_true',
        default=False)
    graph.add_argument(
        "--no_clean_edges",
        dest="clean_edges",
        help=("Turn off edge filtering in the final output graph."),
        action='store_false',
        default=True)

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

    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--quiet",
                        dest="verbose",
                        help="suppress additional output",
                        action='store_false',
                        default=True)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args(args)
    args = set_default_args(args)
    return (args)


def main():
    args = get_options(sys.argv[1:])
    # Check cd-hit is installed
    check_cdhit_version()
    #Make sure aligner is installed if alignment requested
    if args.aln != None:
        check_aligner_install(args.alr)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # check if input is a file containing filenames
    if len(args.input_files) == 1:
        files = []
        with open(args.input_files[0], 'r') as infile:
            for line in infile:
                files.append(line.strip())
        args.input_files = files

    if args.verbose:
        print("pre-processing gff3 files...")

    # convert input GFF3 files into summary files
    process_prokka_input(args.input_files, args.output_dir, (not args.verbose),
                         args.n_cpu)

    # Cluster protein sequences using cdhit
    cd_hit_out = args.output_dir + "combined_protein_cdhit_out.txt"
    run_cdhit(input_file=args.output_dir + "combined_protein_CDS.fasta",
              output_file=cd_hit_out,
              id=args.id,
              s=args.len_dif_percent,
              quiet=(not args.verbose),
              n_cpu=args.n_cpu)

    if args.verbose:
        print("generating initial network...")

    # generate network from clusters and adjacency information
    G, centroid_contexts, seqid_to_centroid = generate_network(
        cluster_file=cd_hit_out + ".clstr",
        data_file=args.output_dir + "gene_data.csv",
        prot_seq_file=args.output_dir + "combined_protein_CDS.fasta",
        all_dna=args.all_seq_in_graph)

    # merge paralogs
    if args.verbose:
        print("Processing paralogs...")
    G = collapse_paralogs(G, centroid_contexts, quiet=(not args.verbose))

    # write out pre-filter graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['genomeIDs'] = ";".join(G.nodes[node]['members'])
        G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
        G.nodes[node]['degrees'] = G.degree[node]
    for edge in G.edges():
        G.edges[edge[0],
                edge[1]]['genomeIDs'] = ";".join(G.edges[edge[0],
                                                         edge[1]]['members'])
    nx.write_gml(G,
                 args.output_dir + "pre_filt_graph.gml",
                 stringizer=custom_stringizer)

    if args.verbose:
        print("collapse mistranslations...")

    # clean up translation errors
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          dna_error_threshold=0.98,
                          correct_mistranslations=True,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose))[0]

    if args.verbose:
        print("collapse gene families...")

    # collapse gene families
    G, distances_bwtn_centroids, centroid_to_index = collapse_families(
        G,
        seqid_to_centroid=seqid_to_centroid,
        outdir=temp_dir,
        family_threshold=args.family_threshold,
        correct_mistranslations=False,
        n_cpu=args.n_cpu,
        quiet=(not args.verbose))

    if args.verbose:
        print("trimming contig ends...")

    # re-trim low support trailing ends
    G = trim_low_support_trailing_ends(G,
                                       min_support=args.min_trailing_support,
                                       max_recursive=args.trailing_recursive)

    if args.verbose:
        print("refinding genes...")

    # find genes that Prokka has missed
    G = find_missing(G,
                     args.input_files,
                     dna_seq_file=args.output_dir + "combined_DNA_CDS.fasta",
                     prot_seq_file=args.output_dir +
                     "combined_protein_CDS.fasta",
                     gene_data_file=args.output_dir + "gene_data.csv",
                     remove_by_consensus=args.remove_by_consensus,
                     search_radius=args.search_radius,
                     prop_match=args.refind_prop_match,
                     pairwise_id_thresh=args.id,
                     merge_id_thresh=max(0.8, args.family_threshold),
                     n_cpu=args.n_cpu,
                     verbose=args.verbose)

    # remove edges that are likely due to misassemblies (by consensus)

    # merge again in case refinding has resolved issues
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          family_threshold=args.family_threshold,
                          correct_mistranslations=False,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose),
                          distances_bwtn_centroids=distances_bwtn_centroids,
                          centroid_to_index=centroid_to_index)[0]

    if args.clean_edges:
        G = clean_misassembly_edges(
            G, edge_support_threshold=args.edge_support_threshold)

    # if requested merge paralogs
    if args.merge_paralogs:
        G = merge_paralogs(G)

    isolate_names = [
        os.path.splitext(os.path.basename(x))[0] for x in args.input_files
    ]
    G.graph['isolateNames'] = isolate_names
    mems_to_isolates = {}
    for i, iso in enumerate(isolate_names):
        mems_to_isolates[str(i)] = iso

    if args.verbose:
        print("writing output...")

    # write out roary like gene_presence_absence.csv
    # get original annotaiton IDs, lengts and whether or
    # not an internal stop codon is present
    orig_ids = {}
    ids_len_stop = {}
    with open(args.output_dir + "gene_data.csv", 'r') as infile:
        next(infile)
        for line in infile:
            line = line.split(",")
            orig_ids[line[2]] = line[3]
            ids_len_stop[line[2]] = (len(line[4]), "*" in line[4][1:-3])

    G = generate_roary_gene_presence_absence(G,
                                             mems_to_isolates=mems_to_isolates,
                                             orig_ids=orig_ids,
                                             ids_len_stop=ids_len_stop,
                                             output_dir=args.output_dir)
    #Write out presence_absence summary
    generate_summary_stats(output_dir=args.output_dir)

    # write pan genome reference fasta file
    generate_pan_genome_reference(G,
                                  output_dir=args.output_dir,
                                  split_paralogs=False)

    # write out common structural differences in a matrix format
    generate_common_struct_presence_absence(
        G,
        output_dir=args.output_dir,
        mems_to_isolates=mems_to_isolates,
        min_variant_support=args.min_edge_support_sv)

    # add helpful attributes and write out graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['centroid'] = ";".join(G.nodes[node]['centroid'])
        G.nodes[node]['dna'] = ";".join(conv_list(G.nodes[node]['dna']))
        G.nodes[node]['protein'] = ";".join(conv_list(
            G.nodes[node]['protein']))
        G.nodes[node]['genomeIDs'] = ";".join(G.nodes[node]['members'])
        G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
        G.nodes[node]['degrees'] = G.degree[node]
        G.nodes[node]['members'] = list(G.nodes[node]['members'])
        G.nodes[node]['seqIDs'] = list(G.nodes[node]['seqIDs'])

    for edge in G.edges():
        G.edges[edge[0],
                edge[1]]['genomeIDs'] = ";".join(G.edges[edge[0],
                                                         edge[1]]['members'])
        G.edges[edge[0],
                edge[1]]['members'] = list(G.edges[edge[0],
                                                   edge[1]]['members'])

    nx.write_gml(G, args.output_dir + "final_graph.gml")

    #Write out core/pan-genome alignments
    if args.aln == "pan":
        if args.verbose: print("generating pan genome MSAs...")
        generate_pan_genome_alignment(G, temp_dir, args.output_dir, args.n_cpu,
                                      args.alr, isolate_names)
        core_nodes = get_core_gene_nodes(G, args.core, len(args.input_files))
        concatenate_core_genome_alignments(core_nodes, args.output_dir)
    elif args.aln == "core":
        if args.verbose: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, args.output_dir,
                                       args.n_cpu, args.alr, isolate_names,
                                       args.core, len(args.input_files))

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
