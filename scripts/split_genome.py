import sys, os
import argparse
from mungo.fasta import FastaReader
import numpy as np


def split_genome(fastafile, outputprefix, n_splits, n_reps):
    seqs = []
    total_length = 0
    for h, s in FastaReader(fastafile):
        seqs.append((h, s))

    if len(seqs) > 1:
        raise ValueError(
            'Multiple contigs in fasta file which isnt currently supported!')
    h = ""
    s = seqs[0][1]

    print "Genome length:", len(s)

    for r in range(n_reps):
        outputfile = outputprefix + "_rep_" + str(r) + ".fasta"
        # generate random split locations
        cuts = np.random.choice(len(s), n_splits)
        cuts = sorted(np.append(cuts, len(s)))
        print "cuts:", cuts
        prev = 0
        with open(outputfile, 'w') as outfile:
            for i, cut in enumerate(cuts):
                outfile.write(">split_" + str(i) + "\n")
                outfile.write(s[prev:cut] + "\n")
                prev = cut

    print "Final contig end position:", len(s)

    return


def main():

    parser = argparse.ArgumentParser(
        description='Fragments a reference genome.')

    parser.add_argument('-i',
                        '--input',
                        dest='input',
                        type=str,
                        required=True,
                        help='input fasta file name')

    parser.add_argument('-o',
                        '--out',
                        dest='out',
                        type=str,
                        required=True,
                        help='output file name (without extension)')

    parser.add_argument(
        '-c',
        '--count',
        dest='n_splits',
        required=True,
        type=int,
        help='number of splits to make. Will result in n_splits + 1 fragments.'
    )

    parser.add_argument('-r',
                        '--rep',
                        dest='n_reps',
                        required=True,
                        type=int,
                        help='number of replicates to run.')

    args = parser.parse_args()

    split_genome(args.input, args.out, args.n_splits, args.n_reps)

    return


if __name__ == '__main__':
    main()
