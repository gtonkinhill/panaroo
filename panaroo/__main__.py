import os, sys
import tempfile
from Bio import SeqIO
import shutil
import networkx as nx
import argparse
import textwrap
import ast

from .isvalid import *
from .set_default_args import set_default_args
from .prokka import process_prokka_input, create_temp_gff3
from .cdhit import check_cdhit_version
from .cdhit import run_cdhit
from .generate_network import generate_network
from .generate_output import *
from .clean_network import *
from .find_missing import find_missing
from .generate_alignments import check_aligner_install
from intbitset import intbitset

from .__init__ import __version__


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            lines = []
            for l in text[2:].splitlines():
                if l == "":
                    lines += [""]
                else:
                    lines += textwrap.wrap(l, width=55)
            return lines
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def get_options(args):

    description = 'panaroo: an updated pipeline for pangenome investigation'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo',
                                     formatter_class=SmartFormatter)

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
                         type=str)

    mode_opts = parser.add_argument_group('Mode')

    mode_opts.add_argument(
        "--clean-mode",
        dest="mode",
        help=
        ('''R|The stringency mode at which to run panaroo. Must be one of 'strict',\
'moderate' or 'sensitive'. Each of these modes can be fine tuned using the\
 additional parameters in the 'Graph correction' section.

strict: 
Requires fairly strong evidence (present in  at least 5%% of genomes)\
 to keep likely contaminant genes. Will remove genes that are refound more often than\
 they were called originally.

moderate: 
Requires moderate evidence (present in  at least 1%% of genomes)\
 to keep likely contaminant genes. Keeps genes that are refound more often than\
 they were called originally.

sensitive: 
Does not delete any genes and only performes merge and refinding\
 operations. Useful if rare plasmids are of interest as these are often hard to\
 disguish from contamination. Results will likely include  higher number of\
 spurious annotations.'''),
        choices=['strict', 'moderate', 'sensitive'],
        required=True)

    mode_opts.add_argument(
        "--remove-invalid-genes",
        dest="filter_invalid",
        action='store_true',
        default=False,
        help=(
            "removes annotations that do not conform to the expected Prokka" +
            " format such as those including premature stop codons."))

    matching = parser.add_argument_group('Matching')
    matching.add_argument("-c",
                          "--threshold",
                          dest="id",
                          help="sequence identity threshold (default=0.98)",
                          type=float)
    matching.add_argument(
        "-f",
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        type=float)
    matching.add_argument("--len_dif_percent",
                          dest="len_dif_percent",
                          help="length difference cutoff (default=0.98)",
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
            "minimum support required to keep an edge that has been flagged" +
            " as a possible mis-assembly"),
        type=float)
    graph.add_argument(
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
    core.add_argument(
        "--codons",
        dest="codons",
        help=
        "Generate codon alignments by aligning sequences at the protein level",
        action='store_true',
        default=False)
    core.add_argument("--core_threshold",
                      dest="core",
                      help="Core-genome sample threshold (default=0.95)",
                      type=float,
                      default=0.95)
    core.add_argument("--core_entropy_filter",
                      dest="hc_threshold",
                      help=("Manually set the Block Mapping and Gathering with " +
                            "Entropy (BMGE) filter. Can be between 0.0 and 1.0. By " + 
                            "default this is set using the Tukey outlier method."),
                      type=float,
                      default=None)

    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--codon-table",
                        dest="table",
                        help="the codon table to use for translation (default=11)",
                        type=int,
                        default=11)
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

    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # check if input is a file containing filenames
    if len(args.input_files) == 1:
        files = []
        with open(args.input_files[0], 'r') as infile:
            for line in infile:
                line = line.strip().split()
                if len(line)==1:
                    ext = os.path.splitext(line[0])[1]
                    print(ext)
                    if ext in ['.gbk', '.gb', '.gbff']:
                        files.append(create_temp_gff3(line[0], None, temp_dir)) 
                    else:  
                        files.append(line[0])
                elif len(line)==2:
                    files.append(create_temp_gff3(line[0], line[1], temp_dir))
                else:
                    print("Problem reading input line: ", line)
                    raise RuntimeError("Error reading files!")
        args.input_files = files

    if args.verbose:
        print("pre-processing gff3 files...")

    # convert input GFF3 files into summary files
    process_prokka_input(args.input_files, args.output_dir,
                         args.filter_invalid, (not args.verbose), 
                         args.n_cpu, args.table)

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
        G.nodes[node]['genomeIDs'] = ";".join(
            [str(m) for m in G.nodes[node]['members']])
        G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
        G.nodes[node]['degrees'] = G.degree[node]
    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            [str(m) for m in G.edges[edge[0], edge[1]]['members']])
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
                          length_outlier_support_proportion=args.
                          length_outlier_support_proportion,
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
        length_outlier_support_proportion=args.
        length_outlier_support_proportion,
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
    if args.verbose:
        print("collapse gene families with refound genes...")
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          family_threshold=args.family_threshold,
                          correct_mistranslations=False,
                          length_outlier_support_proportion=args.
                          length_outlier_support_proportion,
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
        mems_to_isolates[i] = iso

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

    nx.write_gml(G, args.output_dir + "final_graph.gml")

    #Write out core/pan-genome alignments
    if args.aln == "pan":
        if args.verbose: print("generating pan genome MSAs...")
        generate_pan_genome_alignment(G, temp_dir, args.output_dir, args.n_cpu,
                                      args.alr, args.codons, isolate_names)
        core_nodes = get_core_gene_nodes(G, args.core, len(args.input_files))
        core_names = [G.nodes[x]["name"] for x in core_nodes]
        concatenate_core_genome_alignments(core_names, args.output_dir, args.hc_threshold)
    elif args.aln == "core":
        if args.verbose: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, args.output_dir,
                                       args.n_cpu, args.alr, isolate_names,
                                       args.core, args.codons, len(args.input_files),
                                       args.hc_threshold)

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
