#!/usr/bin/env python3
"""
A script to merge a single GFF file with a pre-existing panaroo output
"""

import os
import tempfile
import networkx as nx
import shutil
import sys
import subprocess

from .__init__ import __version__
from .prokka import process_prokka_input
from .cdhit import run_cdhit
from .generate_network import generate_network
from .isvalid import *
from .merge_graphs import merge_graphs


def get_options(
):  #options for integrating (combination of merge graph and cdhit options

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

    io_opts.add_argument("-i",
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
                         type=str)
    
    io_opts.add_argument(
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
    core.add_argument("--core_subset",
                      dest="subset",
                      help="Randomly subset the core genome to these many genes (default=all)",
                      type=int,
                      default=None)
    core.add_argument("--core_entropy_filter",
                      dest="hc_threshold",
                      help=("Manually set the Block Mapping and Gathering with " +
                            "Entropy (BMGE) filter. Can be between 0.0 and 1.0. By " + 
                            "default this is set using the Tukey outlier method."),
                      type=float,
                      default=None)

    graph = parser.add_argument_group('Graph correction')

    graph.add_argument(
        "--all_seq_in_graph",
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
    parser.add_argument("--codon-table",
                        dest="table",
                        help="the codon table to use for translation (default=11)",
                        type=int,
                        default=11)
    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False)

    parser.add_argument(
        "--dirty",
        dest="dirty",
        help=
        "keep temporary directory containing cluster files and cdhit output",
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


def reformat_network(single_gml, output_dir, isolateName):
    """Reformats the output of generate_network() for linear graphs to allow input into merge_graphs()"""
    for adj in single_gml._adj:
        for x in single_gml._adj[adj]:
            y = single_gml._adj[adj][x]

            y.pop('members')
            zero = {'members': 0}
            y.update(zero)

            genomes = {'genomeIDs': '0'}
            y.update(genomes)

    for node in single_gml._node:
        y = single_gml._node[node]
        y.pop('members')

        zero = {
            'members': 0
        }  #members are assigned intbitset[0]. needs to be 0
        y.update(zero)

        to_replace = {"[": "", "]": "", "'": ""}
        y['centroid'] = replace_all(str(y['centroid']), to_replace)
        y['dna'] = replace_all(str(y['dna']), to_replace)
        y['protein'] = replace_all(str(y['protein']), to_replace)

        y['hasEnd'] = int(y['hasEnd'])
        y['mergedDNA'] = int(y['mergedDNA'])
        y['paralog'] = int(y['paralog'])

        y['longCentroidID'] = list(y['longCentroidID'])
        y['seqIDs'] = list(y['seqIDs'])

    single_gml.graph.update({'isolateNames':
                             'x'})  # isolateName from gff filename
    nx.write_gml(single_gml, output_dir + "final_graph.gml")

    return single_gml


def main():
    #Takes a single GFF input, generates a graph and merges with a pre-existing graph
    args = get_options()

    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    args.input_dir = os.path.join(args.input_dir, "")
    args.output_dir = os.path.join(args.output_dir, "")

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    os.environ['TMPDIR'] = temp_dir

    directories = [args.input_dir, temp_dir]

    gff_file = [args.input_gff]

    filename = os.path.basename(args.input_gff).split(".")[0]

    if not args.quiet: print("Processing input")
    process_prokka_input(gff_list=gff_file,
                         output_dir=temp_dir,
                         filter_seqs=args.filter_invalid,
                         quiet=args.quiet,
                         n_cpu=args.n_cpu,
                         table=args.table)

    cd_hit_out = temp_dir + "combined_protein_cdhit_out.txt"

    run_cdhit(input_file=temp_dir + "combined_protein_CDS.fasta",
              output_file=cd_hit_out,
              id=args.id,
              quiet=args.quiet,
              n_cpu=args.n_cpu)

    if not args.quiet: print("Generating network")
    single_gml, centroid_contexts_single, seqid_to_centroid_single = generate_network(
        cluster_file=cd_hit_out + ".clstr",
        data_file=temp_dir + "gene_data.csv",
        prot_seq_file=temp_dir + "combined_protein_CDS.fasta",
        all_dna=args.all_seq_in_graph)

    if not args.quiet: print("Reformatting network")
    reformat_network(single_gml=single_gml,
                     output_dir=temp_dir,
                     isolateName=filename)

    merge_graphs(directories=directories,
                 temp_dir=temp_dir,
                 len_dif_percent=args.len_dif_percent,
                 pid=args.id,
                 family_threshold=args.family_threshold,
                 length_outlier_support_proportion=args.
                 length_outlier_support_proportion,
                 merge_para=args.merge_paralogs,
                 output_dir=args.output_dir,
                 min_edge_support_sv=args.min_edge_support_sv,
                 aln=args.aln,
                 alr=args.alr,
                 core=args.core,
                 codons=args.codons,
                 hc_threshold=args.hc_threshold,
                 subset=args.subset,
                 merge_single=True,
                 depths=[1],
                 n_cpu=args.n_cpu,
                 quiet=args.quiet)

    G = nx.read_gml(args.output_dir + "final_graph.gml")

    for index, name in enumerate(
            G.graph['isolateNames']
    ):  #Corrects isolate name for single gff being returned as list
        if name == 'x':
            G.graph['isolateNames'][index] = filename

    nx.write_gml(G, args.output_dir + "final_graph.gml")

    #remove temporary directory if dirty = True
    if not args.dirty:
        shutil.rmtree(temp_dir)

    sys.exit(0)


if __name__ == '__main__':
    main()
