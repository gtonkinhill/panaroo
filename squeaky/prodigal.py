import subprocess

def run_prodigal(trans_file, nuc_file, output_type, input_file, output_file,
    mode="single",
    shine_dalgarno=True,
    mask_sequence=False,
    tr_table=11,
    closed_ends=False,
    output_type="gff",
    quiet=False,
    start_file=None,
    training_file=None):

    cmd = ("prodigal" +
        " -a " + trans_file +
        " -d " + nuc_file +
        " -f " + output_type +
        " -g " + str(tr_table) +
        " -i " + input_file +
        " -o " + output_file +
        " -p " + mode)

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
