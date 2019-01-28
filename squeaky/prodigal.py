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



# prodigal.run_prodigal(trans_file="./test_processed_files/test_trans_file.txt", nuc_file="./test_processed_files/test_nuc_file.txt", input_file="../test_data/mycobacterium_tuberculosis_H37Rv/GCA_000195955.2_ASM19595v2_genomic.fna", output_file="./test_processed_files/test_output_file.txt")
