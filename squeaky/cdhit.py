import subprocess

def run_cdhit(input_file, output_file,
    id=0.95,
    n_cpu=1,
    aL=0.0, # alignment coverage for the longer sequence
    AL=99999999, # alignment coverage control for the longer sequence
    aS=0.0, # alignment coverage for the shorter sequence
    AS=99999999, # alignment coverage control for the shorter sequence
    accurate=True, # use the slower but more accurate options
    use_local=False, #whether to use local or global sequence alignment
    quiet=False):

    cmd = ("cd-hit" +
        " -T " + str(n_cpu) +
        " -i " + input_file +
        " -o " + output_file +
        " -c " + str(id) +
        " -aL " + str(aL) +
        " -AL " + str(AL) +
        " -aS " + str(aS) +
        " -AS " + str(AS) +
        " -M 0 -d 999"
        )

    if use_local:
        cmd += " -G 0"

    if accurate:
        cmd += " -g 1 -n 2"

    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    return
