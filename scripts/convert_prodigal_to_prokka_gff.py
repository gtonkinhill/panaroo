import sys, os
import argparse
from mungo.fasta import FastaReader
import numpy as np


def convert(gfffile, fastafile, outputfile):

    with open(outputfile, 'w') as outfile:
        outfile.write("##gff-version 3\n")
        for h,s in FastaReader(fastafile):
            h=h.split()[0]
            outfile.write(" ".join(["##sequence-region", h, "1", str(len(s))]) + "\n")

        with open(gfffile, 'rU') as infile:
            for line in infile:
                if line[0]!="#":
                    outfile.write(line)

        outfile.write("##FASTA\n")

        for h,s in FastaReader(fastafile):
            h=h.split()[0]
            outfile.write(">"+h+"\n"+s+"\n")


    return

def main():

    parser = argparse.ArgumentParser(description='Fragments a reference genome.')

    parser.add_argument('-g', '--gff', dest='gff', type=str, required=True,
                       help='input gff file name')

    parser.add_argument('-f', '--fasta', dest='fasta', type=str, required=True,
                       help='input fasta file name')

    parser.add_argument('-o', '--out', dest='out', type=str, required=True,
                       help='output file name')


    args = parser.parse_args()

    convert(args.gff, args.fasta, args.out)

    return



if __name__ == '__main__':
    main()
