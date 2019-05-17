from panaurus.prokka import process_prokka_input
from panaurus.cdhit import run_cdhit
from panaurus.generate_network import generate_network
from panaurus.generate_output import *
from panaurus.clean_network import *
from panaurus.find_missing import find_missing
import os
import argparse
import tempfile
from Bio import SeqIO
import shutil
import networkx as nx
from panaurus.isvalid import *
from panaurus.set_default_args import set_default_args


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-c",
        "--threshold",
        dest="id",
        help="sequence identity threshold (default=0.95)",
        type=float)

    parser.add_argument(
        "-f",
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        type=float)

    parser.add_argument(
        "--len_dif_percent",
        dest="len_dif_percent",
        help="length difference cutoff (default=0.95)",
        type=float)

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="input GFF3 files (usually output from running Prokka)",
        type=argparse.FileType('rU'),
        nargs='+')

    parser.add_argument(
        "-o",
        "--out_dir",
        dest="output_dir",
        required=True,
        help="location of an output directory",
        type=lambda x: is_valid_folder(parser, x))

    parser.add_argument(
        "--min_trailing_support",
        dest="min_trailing_support",
        help=("minimum cluster size to keep a gene called at the " +
              "end of a contig (default=2)"),
        type=int)

    parser.add_argument(
        "--trailing_recursive",
        dest="trailing_recursive",
        help=("number of times to perform recursive triming of low support " +
              "nodes near the end of contigs (default=2)"),
        type=int)

    parser.add_argument(
        "--max_cycle_size",
        dest="max_cycle_size",
        help=("maximum cycle  size for collapsing gene families " +
              "(default=20)"),
        type=int)

    parser.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=("minimum edge support required to call structural variants" +
              " in the presence/absence sv file"),
        type=int)

    parser.add_argument(
        "--no_split",
        dest="split_paralogs",
        help="don't split paralogs",
        action='store_false',
        default=True)

    parser.add_argument(
        "--mode",
        dest="mode",
        help=("the stringency mode at which to run panaurus. One of 'strict'" +
              ", 'moderate' or 'relaxed' (default='strict')"),
        choices=['strict', 'moderate', 'relaxed'],
        default='strict')

    parser.add_argument(
        "-t",
        "--threads",
        dest="n_cpu",
        help="number of threads to use (default=1)",
        type=int,
        default=1)

    parser.add_argument(
        "-a",
        "--alignment",
        dest="aln",
        help=("Output alignments of core genes or all genes. Options are" +
            " 'core' and 'pan'. Default: 'None'"),
        type=str,
        default=None)

    parser.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'prank', 'clustal', and default: 'mafft'",
        type=str,
        default="mafft")

    parser.add_argument(
        "--core_threshold",
        dest="core",
        help="Core-genome sample threshold (default=0.95)",
        type=float,
        default=0.95)

    parser.add_argument(
        "--verbose",
        dest="verbose",
        help="print additional output",
        action='store_true',
        default=False)

    args = parser.parse_args()

    args = set_default_args(args)

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # convert input GFF3 files into summary files
    process_prokka_input(args.input_files, args.output_dir)

    # Cluster protein sequences using cdhit
    cd_hit_out = args.output_dir + "combined_protein_cdhit_out.txt"
    run_cdhit(
        input_file=args.output_dir + "combined_protein_CDS.fasta",
        output_file=cd_hit_out,
        id=args.id,
        s=args.len_dif_percent,
        n_cpu=args.n_cpu)

    if args.verbose:
        print("generating initial network...")

    # generate network from clusters and adjacency information
    G = generate_network(
        cluster_file=cd_hit_out + ".clstr",
        data_file=args.output_dir + "gene_data.csv",
        prot_seq_file=args.output_dir + "combined_protein_CDS.fasta",
        split_paralogs=args.split_paralogs)

    # write out pre-filter graph in GML format
    nx.write_gml(G, args.output_dir + "pre_filt_graph.gml")

    if args.verbose:
        print("triming contig ends...")

    # remove low support trailing ends
    G = trim_low_support_trailing_ends(
        G,
        min_support=args.min_trailing_support,
        max_recursive=args.trailing_recursive)

    if args.verbose:
        print("collapse gene families...")

    # clean up translation errors and gene families
    G = collapse_families(
        G,
        cycle_threshold=args.max_cycle_size,
        family_threshold=args.family_threshold,
        outdir=temp_dir,
        dna_error_threshold=0.99,
        correct_mistranslations=True)

    if args.verbose:
        print("refinding genes...")

    # find genes that Prokka has missed
    G = find_missing(
        G,
        args.input_files,
        temp_dir=temp_dir,
        dna_seq_file=args.output_dir + "combined_DNA_CDS.fasta",
        prot_seq_file=args.output_dir + "combined_protein_CDS.fasta",
        n_cpu=args.n_cpu)

    # write out roary like gene_presence_absence.csv
    G = generate_roary_gene_presence_absence(
        G,
        file_names=args.input_files,
        dna_file=args.output_dir + "combined_DNA_CDS.fasta",
        output_dir=args.output_dir)

    # add helpful attributes and write out graph in GML format
    for node in G.nodes():
        G.node[node]['genomeIDs'] = ";".join(G.node[node]['members'])
        G.node[node]['geneIDs'] = ";".join(G.node[node]['seqIDs'])
        G.node[node]['degrees'] = G.degree[node]
    nx.write_gml(G, args.output_dir + "final_graph.gml")

    # write pan genome reference fasta file
    generate_pan_genome_reference(
        G, output_dir=args.output_dir, split_paralogs=False)

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
        generate_pan_genome_alignment(G, temp_dir, args.output_dir, args.n_cpu, args.alr,
                                      isolate_names)
        core_nodes = get_core_gene_nodes(G, args.core, len(args.input_files))
        concatenate_core_genome_alignments(core_nodes,
                                           args.output_dir)
    elif args.aln == "core":
        if args.verbose: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, args.output_dir, args.n_cpu, args.alr,
                                       isolate_names, args.core,
                                       len(args.input_files))

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
