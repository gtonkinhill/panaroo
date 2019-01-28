import subprocess

def run_jackhmmer(seqfile, seqdb, output_file,
    n_iter=5,
    n_cpu=1,
    E=None,
    T=None,
    domE=None,
    domT=None,
    incE=None,
    incT=None,
    incdomE=None,
    incdomT=None,
    align_file=None,
    tblout=None,
    domtblout=None,
    quiet=False):

    cmd = ("jackhmmer" +
        " --cpu " + n_cpu +
        " -N " + n_iter +
        " -o " + output_file)

    if not E is None:
        cmd += " -E " + E

    if not T is None:
        cmd += " -T " + T

    if not domE is None:
        cmd += " --domE " + domE

    if not domT is None:
        cmd += " --domT " + domT

    if not incE is None:
        cmd += " --incE " + incE

    if not incT is None:
        cmd += " --incT " + incT

    if not incdomE is None:
        cmd += " --incdomE " + incdomE

    if not incdomT is None:
        cmd += " --incdomT " + incdomT

    if not align_file is None:
        cmd += " -A " + align_file

    if not tblout is None:
        cmd += " --tblout " + tblout

    if not domtblout is None:
        cmd += " --domtblout " + domtblout

    cmd += seqfile + " " + seqdb

    if not quiet:
        print("running cmd: " + cmd)


    subprocess.run(cmd, shell=True, check=True)

    return
