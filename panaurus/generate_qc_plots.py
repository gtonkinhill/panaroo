import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import subprocess
import tempfile
import os

from .isvalid import *
from .__init__ import __version__


def generate_plot(input_gffs, method, outdir):

    return


def get_options():
    import argparse

    description = 'Generate summary plots from Panaurus run'
    parser = argparse.ArgumentParser(description=description,
                                     prog='plot_panaurus_summary')

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

    # Other options
    parser.add_argument("--graph_type",
                        dest="graph_type",
                        help="the type of graph to generate",
                        choices={'all'})
    parser.add_argument("--ref_db",
                        dest="ref_db",
                        help="reference mash database for contamination quantification.")
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    if args.graph_type == 'all':
        for method in []:
            generate_plot(args.input_files, method, args.output_dir)
    else:
        generate_plot(args.input_files, args.graph_type, args.output_dir)

    return


if __name__ == '__main__':
    main()
