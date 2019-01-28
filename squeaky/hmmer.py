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
        " --cpu " + str(n_cpu) +
        " -N " + str(n_iter) +
        " -o " + output_file)

    if not E is None:
        cmd += " -E " + str(E)

    if not T is None:
        cmd += " -T " + str(T)

    if not domE is None:
        cmd += " --domE " + str(domE)

    if not domT is None:
        cmd += " --domT " + str(domT)

    if not incE is None:
        cmd += " --incE " + str(incE)

    if not incT is None:
        cmd += " --incT " + str(incT)

    if not incdomE is None:
        cmd += " --incdomE " + str(incdomE)

    if not incdomT is None:
        cmd += " --incdomT " + str(incdomT)

    if not align_file is None:
        cmd += " -A " + align_file

    if not tblout is None:
        cmd += " --tblout " + tblout

    if not domtblout is None:
        cmd += " --domtblout " + domtblout

    cmd += " " + seqfile + " " + seqdb

    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    return


# hmmer.run_jackhmmer("./test_processed_files/test_query.fasta", "./test_processed_files/test_trans_file.fasta", "./test_processed_files/test_hmmer.txt")
