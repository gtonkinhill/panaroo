from prodigal import run_prodigal
from shutil import copyfileobj
import os
import tempfile
from joblib import Parallel, delayed
from Bio import SeqIO


def run_prodigal_combined(input_file, training_file, temp_dir,
    mode="single",
    shine_dalgarno=True,
    mask_sequence=False,
    tr_table=11,
    closed_ends=False,
    output_type="gff",
    quiet=False,
    start_file=None):

    temp_trans = tempfile.NamedTemporaryFile(delete = False, dir=temp_dir)
    temp_nuc = tempfile.NamedTemporaryFile(delete = False, dir=temp_dir)
    temp_out = tempfile.NamedTemporaryFile(delete = False, dir=temp_dir)
    temp_trans.close()
    temp_nuc.close()
    temp_out.close()

    run_prodigal(trans_file=temp_trans.name, nuc_file=temp_nuc.name,
        input_file=input_file, output_file=temp_out.name,
        training_file=training_file,
        mode=mode, shine_dalgarno=shine_dalgarno,
        mask_sequence=mask_sequence, tr_table=tr_table,
        closed_ends=closed_ends,
        quiet=quiet)

    return (temp_trans.name, temp_nuc.name, temp_out.name)


def run_prodigal_multi(input_files, out_dir, temp_dir,
    mode="single",
    shine_dalgarno=True,
    mask_sequence=False,
    tr_table=11,
    closed_ends=False,
    output_type="gff",
    quiet=False,
    start_file=None,
    training_file=None,
    n_cpu=1,
    max_n_training=100):


    # take first 100 files for training prodigal
    temp_file = tempfile.NamedTemporaryFile(delete = False, dir=temp_dir)
    for f in input_files[:max_n_training]:
        with open(f, 'rb') as inputfile:
            copyfileobj(inputfile, temp_file)
    temp_file.close()

    # generate training file
    training_file = out_dir + "prodigal_training_file.txt"
    if not os.path.exists(training_file):
        run_prodigal(trans_file="/dev/null", nuc_file="/dev/null",
            input_file=temp_file.name, output_file="/dev/null",
            training_file=training_file,
            mode=mode, shine_dalgarno=shine_dalgarno,
            mask_sequence=mask_sequence, tr_table=tr_table,
            closed_ends=closed_ends,
            quiet=quiet)
    # delete temp file
    os.remove(temp_file.name)

    # run prodigal on each assembly in parallel
    out_file_list = Parallel(n_jobs=n_cpu)(delayed(
        run_prodigal_combined)(f, training_file, temp_dir) for f in input_files)

    # concatenate results into one file and rename
    with open(out_dir + "combined_protein_CDS.fasta", 'a') as outfile:
        for f in out_file_list:
            with open(f[0], 'rU') as inputfile:
                for rec in SeqIO.parse(inputfile, "fasta"):
                    outfile.write(">" + str(rec.id).split()[0] + "\n" +
                        str(rec.seq) + "\n")

    with open(out_dir + "combined_DNA_CDS.fasta", 'a') as outfile:
        for f in out_file_list:
            with open(f[1], 'rU') as inputfile:
                for rec in SeqIO.parse(inputfile, "fasta"):
                    outfile.write(">" + str(rec.id).split()[0] + "\n" +
                        str(rec.seq) + "\n")

    return
