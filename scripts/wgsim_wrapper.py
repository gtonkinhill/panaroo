import subprocess
from joblib import Parallel, delayed


def simulate_reads(reference, read_number, read_length, outname):
    command = "wgsim -N " + str(read_number) + " -1 " + str(
        read_length) + " -2 " + str(read_length) + " "
    command = command + reference + " " + outname + "_1.fastq " + outname + "_2.fastq"
    print(command)
    subprocess.run(command, shell=True)
    return True


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="""Simulate some reads""")
    parser.add_argument('--readLength',
                        help="Provide the length of the read to simulate",
                        default=125)
    parser.add_argument('--threads',
                        help="No. of threads for multithreading, default=1",
                        default=1,
                        type=int)
    parser.add_argument('--readNo',
                        help="Number of reads to simulate",
                        type=int)
    parser.add_argument('--repNo',
                        help="Nubmer of simulations to run",
                        type=int)
    parser.add_argument("reference", help="Reference fasta to simulate reads")

    args = parser.parse_args()
    if args.threads > 1:
        reps = range(args.repNo)
        name = args.reference.split('/')[-1].split('.')[0]
        outnames = [name + "-" + str(x) for x in reps]
        dummy_list = Parallel(n_jobs=args.threads)(delayed(simulate_reads)(
            args.reference, args.readNo, args.readLength, out)
                                                   for out in outnames)
    else:
        print("Single threading not yet supported")
