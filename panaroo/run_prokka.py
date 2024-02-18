import os
import argparse
from .prodigal import train_prodigal
from .isvalid import *
import subprocess
from joblib import Parallel, delayed
import shutil
import tempfile
from tqdm import tqdm

from .__init__ import __version__


def get_options():
    import argparse

    description = 'Re-run gene finding and annotation'
    parser = argparse.ArgumentParser(description=description,
                                     prog='run_prokka')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument("-i",
                         "--input",
                         dest="input_files",
                         required=True,
                         help="input fasta files",
                         type=argparse.FileType('r'),
                         nargs='+')
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=lambda x: is_valid_folder(parser, x))
    io_opts.add_argument("--force",
                         dest="force",
                         help="overwrite old commands",
                         action='store_true',
                         default=False)

    prokka = parser.add_argument_group('Prokka/prodigal command')
    prokka.add_argument(
        "--add_prokka_cmds",
        dest="add_prokka_cmds",
        help="additional commands to supply to Prokka (these are not checked!)",
        type=str,
        default=None)
    prokka.add_argument(
        "--num_training",
        dest="num_training",
        help="number of genomes to use in training prodigal (default=10)",
        type=int,
        default=10)

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

    if len(args.input_files) < args.num_training:
        args.num_training = len(args.input_files)

    # run prodigal to generate training file if it doesn't already exist.
    train_prodigal(input_files=args.input_files,
                   n_samples=args.num_training,
                   force=args.force,
                   outdir=args.output_dir)

    # run prokka with adjusted arguments on each fasta input in parallel
    Parallel(n_jobs=args.n_cpu)(delayed(run_prokka_mod)(
        input, args.output_dir, args.output_dir +
        "prodigal_training.txt", args.force, args.add_prokka_cmds)
                                for input in tqdm(args.input_files))

    return


def run_prokka_mod(input_file, out_folder, train_file, force, add_cmds):

    prefix = os.path.splitext(os.path.basename(input_file.name))[0]
    out_folder = os.path.join(out_folder, "")

    path_to_prodigal = shutil.which("prodigal")

    # override the system prodigal temporarily to input training file in prokka
    # Create temporary directory
    temp_dir = os.path.join(os.path.abspath(tempfile.mkdtemp(dir=out_folder)),
                            "")
    os.environ['TMPDIR'] = temp_dir
    with open(temp_dir + "/prodigal", 'w') as outfile:
        outfile.write("#!/bin/bash\n")
        outfile.write("(>&2 echo 'running prokka mod!')\n")
        outfile.write(path_to_prodigal + " -t " + os.path.abspath(train_file) +
                      " $*\n")

    cmd = 'export PATH="' + temp_dir + ':$PATH"; chmod 777 ' + temp_dir + '/*; '

    cmd += "prokka"

    # TODO: This made be a bad thing to do security wise
    if add_cmds is not None:
        cmd += " " + add_cmds + " "

    if force:
        cmd += " --force"

    cmd += " --cpus 1"
    cmd += " --outdir " + out_folder + prefix
    cmd += " --prefix " + prefix
    cmd += " " + input_file.name
    cmd += "  > " + out_folder + prefix + "_prokka.log 2>&1"

    # print(cmd)

    subprocess.run(cmd, shell=True, check=True, executable='/bin/sh')

    # check prokka completed successfully
    with open(out_folder + prefix + "_prokka.log", 'r') as logfile:
        lines = logfile.read().splitlines()
        if "Annotation finished successfully." not in lines[-6]:
            raise Exception('Prokka did not execute successfully!')

    # remove temporary directory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
