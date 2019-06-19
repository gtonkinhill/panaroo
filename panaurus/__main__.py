import os
import tempfile
from Bio import SeqIO
import shutil
import networkx as nx

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

    description = 'panaurus: an updated pipeline for pan-genome investigation'
    parser = argparse.ArgumentParser(description=description, prog='panaurus')

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

    graph = parser.add_argument_group('Graph correction')
    graph.add_argument(
        "--mode",
        dest="mode",
        help=("the stringency mode at which to run panaurus. One of 'strict'" +
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
        "--edge_support_diff",
        dest="edge_support_diff",
        help=("maximum fraction difference between an edge's support " +
              "and those of the nodes it connects"),
        type=float)
    graph.add_argument(
        "--remove_by_consensus",
        dest="remove_by_consensus",
        help=
        ("if a gene is called in the same region with similar sequence a minority "
         + "of the time, remove it"),
        action='store_true',
        default=None)
    graph.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=("minimum edge support required to call structural variants" +
              " in the presence/absence sv file"),
        type=int)

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
                         "combined_protein_CDS.fasta")

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
        print("collapse gene families...")

    # clean up translation errors
    G = collapse_families(G,
                          outdir=temp_dir,
                          dna_error_threshold=0.95,
                          correct_mistranslations=True,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose))

    # collapse gene families
    G = collapse_families(G,
                          outdir=temp_dir,
                          family_threshold=args.family_threshold,
                          correct_mistranslations=False,
                          n_cpu=args.n_cpu,
                          quiet=(not args.verbose))

    # remove edges that are likely due to misassemblies (by consensus)
    G = clean_misassembly_edges(G, threshold=args.edge_support_diff)

    # re-trim low support trailing ends
    G = trim_low_support_trailing_ends(G,
                                       min_support=args.min_trailing_support,
                                       max_recursive=args.trailing_recursive)

    if args.verbose:
        print("refinding genes...")

    # identify possible family level paralogs
    if args.verbose:
        print("identifying potentialy highly variable genes...")
    G = identify_possible_highly_variable(G)

    # find genes that Prokka has missed
    G = find_missing(G,
                     args.input_files,
                     temp_dir=temp_dir,
                     dna_seq_file=args.output_dir + "combined_DNA_CDS.fasta",
                     prot_seq_file=args.output_dir +
                     "combined_protein_CDS.fasta",
                     remove_by_consensus=args.remove_by_consensus,
                     n_cpu=args.n_cpu)

    # if requested merge paralogs
    if args.merge_paralogs:
        G = merge_paralogs(G)

    # write out roary like gene_presence_absence.csv
    G = generate_roary_gene_presence_absence(G,
                                             file_names=args.input_files,
                                             dna_file=args.output_dir +
                                             "combined_DNA_CDS.fasta",
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
        n_members=len(args.input_files),
        min_variant_support=args.min_edge_support_sv)

    #Write out core/pan-genome alignments
    isolate_names = [
        os.path.splitext(os.path.basename(x.name))[0] for x in args.input_files
    ]
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
