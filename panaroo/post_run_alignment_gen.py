import shutil
import tempfile
import os
import networkx as nx
from .generate_output import *

from .isvalid import *
from .__init__ import __version__


def get_options():
    import argparse

    description = 'Generate multiple sequence alignments after running Panaroo'
    parser = argparse.ArgumentParser(description=description,
                                     prog='generate_panaroo_msa')

    io_opts = parser.add_argument_group('Input/output')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of the Panaroo output directory",
                         type=lambda x: is_valid_folder(parser, x))

    # alignment
    core = parser.add_argument_group('Gene alignment')
    core.add_argument(
        "-a",
        "--alignment",
        dest="aln",
        help=("Output alignments of core genes or all genes. Options are" +
              " 'core' and 'pan'. Default: 'None'"),
        type=str,
        choices={'core', 'pan'},
        default='core')
    core.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'prank', 'clustal', and default: 'mafft'",
        type=str,
        choices={'prank', 'clustal', 'mafft'},
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
    parser.add_argument("--verbose",
                        dest="verbose",
                        help="print additional output",
                        action='store_true',
                        default=False)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # Load isolate names
    seen = set()
    isolate_names = []
    with open(args.output_dir + "gene_data.csv", 'r') as infile:
        next(infile)
        for line in infile:
            iso = line.split(",")[0]
            if iso not in seen:
                isolate_names.append(iso)
                seen.add(iso)

    # Load graph
    G = nx.read_gml(args.output_dir + "final_graph.gml")

    #Write out core/pan-genome alignments
    if args.aln == "pan":
        if args.verbose: print("generating pan genome MSAs...")
        generate_pan_genome_alignment(G, temp_dir, args.output_dir, args.n_cpu,
                                      args.alr, args.codons, isolate_names)

        core_nodes = get_core_gene_nodes(G, args.core, len(isolate_names))
        core_names = [G.nodes[x]["name"] for x in core_nodes]
        concatenate_core_genome_alignments(core_names, args.output_dir)
    elif args.aln == "core":
        if args.verbose: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, args.output_dir,
                                       args.n_cpu, args.alr, isolate_names,
                                       args.core, args.codons, len(isolate_names),
                                       args.hc_threshold)

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
