from prodigal import run_prodigal
from cdhit import run_cdhit
import os
import argparse


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

def is_valid_folder(parser, arg):
    if not os.path.isdir(arg):
        parser.error("The folder %s does not exist!" % arg)
    else:
        return arg



def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--threshold", dest="id",
                    help="sequence identity threshold (default=0.95)",
                    type=float, default=0.95)

    parser.add_argument("-i", "--input", dest="input_file", required=True,
                    help="input file",
                    type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-p", "--train_dir", dest="train_dir",
                    help=("location of a training directory (a previous run" +
                        " of squeaky)"),
                    type=lambda x: is_valid_folder(parser, x))

    parser.add_argument("-o", "--out_dir", dest="output_dir", required=True,
                    help="location of an output directory",
                    type=lambda x: is_valid_folder(parser, x))

    parser.add_argument("-t", "--threads", dest="n_cpu",
                    help="number of threads to use (default=1)",
                    type=int, default=1)

    args = parser.parse_args()

    # Get prefix
    prefix = os.path.splitext(os.path.basename(args.input_file))[0]
    prot_file = args.output_dir + "/" + prefix + "_prot_CDS.fasta"

    # Run prodigal on new sequences
    run_prodigal(trans_file=prot_file,
        nuc_file=args.output_dir + "/" + prefix + "_dna_CDS.fasta",
        input_file=args.input_file,
        output_file=args.output_dir + "/" + prefix + "_prodigal.txt")

    # Cluster protein sequences using cdhit
    run_cdhit(input_file=prot_file,
        output_file=args.output_dir + "/" + prefix + "_cdhit_out.txt",
        id=args.id,
        n_cpu=args.n_cpu)

    return


if __name__ == '__main__':
    main()
