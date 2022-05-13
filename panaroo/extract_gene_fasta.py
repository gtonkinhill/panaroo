# Originally written by Sam Lipworth (https://github.com/samlipworth)
import os, sys
from .isvalid import *
from .__init__ import __version__

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_options(args):
    import argparse

    parser = argparse.ArgumentParser(
        description="Extracts fasta files for a given gene name from a Panaroo pangenome graph.",
        prog="panaroo_extract_fasta",
    )

    parser.add_argument(
        "-q"
        "--queries", 
        dest="queries", 
        required=True, 
        help="A list of gene names to extract",
        type=str,
        nargs='+'
        )

    parser.add_argument(
        "--pa", 
        dest="pa_file", 
        required=True, 
        help="The 'gene_presence_absence.csv' file output by Panaroo",
        type = str
    )

    parser.add_argument(
        "--gene", 
        dest="gene_data", 
        required=True, 
        help="The 'gene_data.csv' file output by Panaroo",
        type = str
    )

    parser.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=lambda x: is_valid_folder(parser, x))


    parser.add_argument("--dna",
                        dest="isdna",
                        help="toggles the output of nucleotide sequence instead of protein sequence",
                        action='store_true',
                        default=False)

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.parse_args(args)
    return (args)


def generate_fasta(geneids, outputfile, genedata, isdna):

    with open(outputfile, 'w') as outfile:
        with open(genedata, 'r') as infile:
            for line in infile:
                line = line.strip().split(',')
                if line[3] in geneids:
                    if isdna:
                        seq = line[5]
                    else:
                        seq = line[4]
                    SeqIO.write(
                        SeqRecord(Seq(seq),
                        id=line[3],
                        description=""), 
                        outfile, "fasta")

    return


def main():
    args = get_options(sys.argv[1:])


    with open(args.pa_file, 'r') as infile:
        header = next(infile).strip().split(',')
        genomes = header[4:]

        for line in infile:
            line = line.strip().split(',')
            if line[0] in args.queries:
                geneids = set([g.replace('_len', '').replace('_stop', '') for g in line[4:]])
                generate_fasta(geneids, 
                    outputfile = args.output_dir + line[0] + '.fasta', 
                    genedata = args.gene_data, 
                    isdna = args.isdna)

    return

if __name__ == "__main__":
    main()