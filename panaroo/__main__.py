import os
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

from .__init__ import __version__


def get_options():
    import argparse

    description = 'panaroo: an updated pipeline for pan-genome investigation'
    parser = argparse.ArgumentParser(description=description, prog='panaroo')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="input GFF3 files (usually output from running Prokka)",
        type=argparse.FileType('rU'),
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
        default=10000,
        type=int)
    refind.add_argument(
        "--refind_prop_match",
        dest="refind_prop_match",
        help=("the proportion of an accessory gene that must " +
              "be found in order to consider it a match"),
        default=0.2,
        type=int)

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
        help=("number of times to perform recursive triming of low support " +
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

    core = parser.add_argument_group('Gene alignment')
    core.add_argument(
        "-a",
        "--alignment",
        dest="aln",
        help=("Output alignments of core genes or all genes. Options are" +
              " 'core' and 'pan'. Default: 'None'"),
        type=str,
        default=None)
    core.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'prank', 'clustal', and default: 'mafft'",
        type=str,
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
    parser.add_argument("--verbose",
                        dest="verbose",
                        help="print additional output",
                        action='store_true',
                        default=False)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    args = set_default_args(args)
    return (args)


def main():
    args = get_options()
    # Check cd-hit is installed
    check_cdhit_version()
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    if args.verbose:
        print("pre-processing gff3 files...")

    # convert input GFF3 files into summary files
    process_prokka_input(args.input_files, args.output_dir, args.n_cpu)

    # Cluster protein sequences using cdhit
    cd_hit_out = args.output_dir + "combined_protein_cdhit_out.txt"
    run_cdhit(input_file=args.output_dir + "combined_protein_CDS.fasta",
              output_file=cd_hit_out,
              id=args.id,
              s=args.len_dif_percent,
              n_cpu=args.n_cpu)

    if args.verbose:
        print("generating initial network...")

    # generate network from clusters and adjacency information
    G = generate_network(cluster_file=cd_hit_out + ".clstr",
                         data_file=args.output_dir + "gene_data.csv",
                         prot_seq_file=args.output_dir +
                         "combined_protein_CDS.fasta",
                         all_dna=args.all_seq_in_graph)

    # merge paralogs
    G = collapse_paralogs(G)

    # write out pre-filter graph in GML format
    for node in G.nodes():
        G.node[node]['size'] = len(set(G.node[node]['members']))
        G.node[node]['genomeIDs'] = ";".join(conv_list(
            G.node[node]['members']))
        G.node[node]['geneIDs'] = ";".join(conv_list(G.node[node]['seqIDs']))
        G.node[node]['degrees'] = G.degree[node]
    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            conv_list(G.edges[edge[0], edge[1]]['members']))
    nx.write_gml(G, args.output_dir + "pre_filt_graph.gml")

    if args.verbose:
        print("triming contig ends...")

    # remove low support trailing ends
    G = trim_low_support_trailing_ends(G,
                                       min_support=args.min_trailing_support,
                                       max_recursive=args.trailing_recursive)

    if args.verbose:
        print("collapse mistranslations...")

    # clean up translation errors
    G = collapse_families(G,
                          outdir=temp_dir,
                          dna_error_threshold=0.99,
                          correct_mistranslations=True,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose))

    if args.verbose:
        print("collapse gene families...")

    # collapse gene families
    G = collapse_families(G,
                          outdir=temp_dir,
                          family_threshold=args.family_threshold,
                          correct_mistranslations=False,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose))

    # re-trim low support trailing ends
    G = trim_low_support_trailing_ends(G,
                                       min_support=args.min_trailing_support,
                                       max_recursive=args.trailing_recursive)

    # identify possible family level paralogs
    if args.verbose:
        print("identifying potentialy highly variable genes...")
    G = identify_possible_highly_variable(
        G,
        cycle_threshold_max=20,
        cycle_threshold_min=args.cycle_threshold_min,
        size_diff_threshold=0.5)

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
                     n_cpu=args.n_cpu)

    # remove edges that are likely due to misassemblies (by consensus)
    G = clean_misassembly_edges(
        G, edge_support_threshold=args.edge_support_threshold)

    # if requested merge paralogs
    if args.merge_paralogs:
        G = merge_paralogs(G)

    isolate_names = [
        os.path.splitext(os.path.basename(x.name))[0] for x in args.input_files
    ]
    G.graph['isolateNames'] = isolate_names
    mems_to_isolates = {}
    for i, iso in enumerate(isolate_names):
        mems_to_isolates[str(i)] = iso

    # write out roary like gene_presence_absence.csv
    G = generate_roary_gene_presence_absence(G,
                                             mems_to_isolates=mems_to_isolates,
                                             output_dir=args.output_dir)

    # add helpful attributes and write out graph in GML format
    for node in G.nodes():
        G.node[node]['size'] = len(set(G.node[node]['members']))

        G.node[node]['genomeIDs'] = ";".join(conv_list(
            G.node[node]['members']))
        G.node[node]['geneIDs'] = ";".join(conv_list(G.node[node]['seqIDs']))
        G.node[node]['degrees'] = G.degree[node]

    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            conv_list(G.edges[edge[0], edge[1]]['members']))

    nx.write_gml(G, args.output_dir + "final_graph.gml")

    # write pan genome reference fasta file
    generate_pan_genome_reference(G,
                                  output_dir=args.output_dir,
                                  split_paralogs=False)

    # write out csv indicating the mobility of each gene
    generate_gene_mobility(G, output_dir=args.output_dir)

    # write out common structural differences in a matrix format
    generate_common_struct_presence_absence(
        G,
        output_dir=args.output_dir,
        mems_to_isolates=mems_to_isolates,
        min_variant_support=args.min_edge_support_sv)

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
