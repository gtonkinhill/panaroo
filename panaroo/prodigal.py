import subprocess
from random import sample
import tempfile
import os


def run_prodigal(trans_file,
                 nuc_file,
                 input_file,
                 output_file,
                 mode="single",
                 shine_dalgarno=True,
                 mask_sequence=False,
                 tr_table=11,
                 closed_ends=False,
                 output_type="gff",
                 quiet=False,
                 start_file=None,
                 training_file=None):
    """Runs prodigal on the command line

    Runs prodigal. Will create a training file if it doesn't already exist.

    Args:
        trans_file (str): output location of translated sequences
        nuc_file (str): output location of nucleotide sequences
        input_file (str): location of input fasta file
        output_file (str): output location of prodigal annotations

    Returns:

    """

    cmd = "prodigal"
    cmd += " -a " + trans_file
    cmd += " -d " + nuc_file
    cmd += " -f " + output_type
    cmd += " -g " + str(tr_table)
    cmd += " -i " + input_file
    cmd += " -o " + output_file
    cmd += " -p " + mode

    if closed_ends:
        cmd += " -c"

    if mask_sequence:
        cmd += " -m"

    if not shine_dalgarno:
        cmd += " -n"

    if not start_file is None:
        cmd += " -s " + start_file

    if not training_file is None:
        cmd += " -t " + training_file

    if quiet:
        cmd += " -q"
    else:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    return


def train_prodigal(input_files, n_samples, force, outdir):

    if os.path.exists(outdir + "prodigal_training.txt"):
        # we already have performed training
        if force:
            os.remove(outdir + "prodigal_training.txt")
        else:
            return

    # randomly sample training genomes
    if len(input_files) < n_samples:
        samples = input_files
    else:
        samples = sample(input_files, n_samples)

    # combined training files into one fasta
    temp_input_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_input_file.close()
    with open(temp_input_file.name, 'w') as outfile:
        for infile in input_files:
            outfile.write(infile.read().strip() + "\n")

    # run prodigal
    run_prodigal(
        trans_file="/dev/null",
        nuc_file="/dev/null",
        input_file=temp_input_file.name,
        output_file="/dev/null",
        training_file=outdir + "prodigal_training.txt",
    )

    os.remove(temp_input_file.name)

    return
