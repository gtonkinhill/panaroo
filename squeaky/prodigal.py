import subprocess

def run_prodigal(trans_file, nuc_file, input_file, output_file,
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

if __name__ == "__main__":
    #For test running prodigal on a single assembly with this script
    from sys import argv
    script, prodigal_input = argv
    filename = prodigal_input.split('/')[-1].split('.')[0]
    run_prodigal(trans_file="./prodigal_test/"+filename+"_protTrans.fasta", nuc_file="./prodigal_test/"+filename+"nucFileTest.fasta", input_file=prodigal_input, output_file="./prodigal_test/"+filename.split(".")[0]+'.gff')
