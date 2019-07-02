import sys, os
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import numpy as np
from random import sample 



def add_contam(reference_fasta, contam_fasta, output_file,
    contam_mean_length, mean_contigs):

    outfile = open(output_file, 'w')

    # write reference to out
    with open(reference_fasta, 'r') as infile:
        outfile.write(infile.read())

    # load contamination
    contam_seqs = list(SeqIO.parse(contam_fasta, 'fasta'))
    contam_seqs = [s for s in contam_seqs if len(s.seq)>500000]

    n_contigs = int(np.random.poisson(lam=mean_contigs, size=1))
    if n_contigs<1: return

    for i in range(n_contigs):
        contam_length = int(np.random.poisson(lam=contam_mean_length, size=1))
        cseq = str(sample(contam_seqs, 1)[0].seq)
        if (len(cseq)<contam_length):
            start=0
        else:
            start = int(np.random.choice(len(cseq)-contam_length, size=1))
        print(start, len(cseq), contam_length)
        cseq = cseq[start:(start+contam_length)]
        outfile.write(">contam_" + str(i) + "\n")
        outfile.write(str(cseq))
    
    # close file
    outfile.close()

    return


def main():

    parser = argparse.ArgumentParser(
        description='Adds artificial variation to genes in a gff.')

    parser.add_argument('-f',
                        '--fasta',
                        dest='fasta',
                        type=str,
                        required=True,
                        help='input reference fasta file name')

    parser.add_argument('-c',
                        '--contam',
                        dest='contam',
                        type=str,
                        required=True,
                        help='input contamination fasta file name')

    parser.add_argument('-l', '--contam_mean_length',
                        dest='contam_mean_length',
                        type=float,
                        required=True,
                        help='mean length of contaminating contigs')

    parser.add_argument('--nreps',
                        dest='nreps',
                        type=int,
                        required=True,
                        help='number of replications')

    parser.add_argument('-m','--mean_contigs',
                        dest='mean_contigs',
                        type=float,
                        required=True,
                        help='mean number of contaminant contigs')

    parser.add_argument('-o',
                        '--out',
                        dest='output_dir',
                        type=str,
                        required=True,
                        help='output directory')

    args = parser.parse_args()
    args.output_dir = os.path.join(args.output_dir, "")

    prefix = os.path.splitext(os.path.basename(args.fasta))[0]
    prefix_contam = os.path.splitext(os.path.basename(args.contam))[0]

    for i in range(args.nreps):
        out_file_name = (args.output_dir + prefix + "_contam_" +
                        prefix_contam + "_rep_" + str(i) + ".fasta")

        add_contam(args.fasta, args.contam, out_file_name, 
            contam_mean_length=args.contam_mean_length,
            mean_contigs=args.mean_contigs)

    return


if __name__ == '__main__':
    main()
