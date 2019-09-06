# A large chunk of this conversion code has been taken from https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py

import argparse
import subprocess
import sys, os
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio import Seq
from io import StringIO

from BCBio import GFF


def _fix_ncbi_id(fasta_iter):
    """GenBank identifiers can only be 16 characters; try to shorten NCBI.
    """
    for rec in fasta_iter:
        if len(rec.name) > 16 and rec.name.find("|") > 0:
            new_id = [x for x in rec.name.split("|") if x][-1]
            print("Warning: shortening NCBI name %s to %s" % (rec.id, new_id))
            rec.id = new_id
            rec.name = new_id
        yield rec

def _check_gff(gff_iterator):
    """Check GFF files before feeding to SeqIO to be sure they have sequences.
    """
    for rec in gff_iterator:
        if isinstance(rec.seq, Seq.UnknownSeq):
            print("Warning: FASTA sequence not found for '%s' in GFF file" % (
                    rec.id))
            rec.seq.alphabet = generic_dna
        yield _flatten_features(rec)

def _flatten_features(rec):
    """Make sub_features in an input rec flat for output.
    GenBank does not handle nested features, so we want to make
    everything top level.
    """
    out = []
    for f in rec.features:
        cur = [f]
        while len(cur) > 0:
            nextf = []
            for curf in cur:
                out.append(curf)
                if len(curf.sub_features) > 0:
                    nextf.extend(curf.sub_features)
            cur = nextf
    rec.features = out
    return rec

#Clean other "##" starting lines from gff file, as it confuses parsers
def clean_gff_string(gff_string):
    splitlines = gff_string.splitlines()
    lines_to_delete = []
    for index in range(len(splitlines)):
        if '##sequence-region' in splitlines[index]:
            lines_to_delete.append(index)
    for index in sorted(lines_to_delete, reverse=True):
        del splitlines[index]
    cleaned_gff = "\n".join(splitlines)
    return cleaned_gff

def convert_to_gbk(gff_file_name, outdir):
    prefix = os.path.splitext(os.path.basename(gff_file_name))[0]
    gff_file = open(gff_file_name, 'r')
    line_list = gff_file.readlines()
    locus_dict = {}
    contig_count = 0
    bad_lines = []

    # deal with long contig names by renaming them
    for i in range(len(line_list)):
        if "##sequence-region" in line_list[i]:    
            locus = line_list[i].split()[1]
            locus_dict[locus] = "Contig" + str(contig_count)
            contig_count+=1
            line_list[i] = line_list[i].replace(locus, locus_dict[locus])
        elif "ID=" in line_list[i]:
            locus = line_list[i].split()[0]
            line_list[i] = line_list[i].replace(locus, locus_dict[locus])
        elif line_list[i][0]==">":
            locus = line_list[i].strip()[1:]
            line_list[i] = line_list[i].replace(locus, locus_dict[locus])
        elif (len(line_list[i].split())>3) and ("=" in line_list[i]):
            locus = line_list[i].split()[0]
            if locus not in locus_dict:
                bad_lines.append(i)
            else:
                line_list[i] = line_list[i].replace(locus, locus_dict[locus])

    for i in bad_lines:
        del line_list[i] 

    #Split file and parse
    lines = "".join(line_list)
    split = lines.split('##FASTA')
    if len(split) != 2:
        print("Problem reading GFF3 file: ", gff_file.name)
        raise RuntimeError("Error reading prokka input!")
    out_file = outdir+prefix+".gbk"
    with StringIO(split[1]) as temp_fasta:
        with StringIO(clean_gff_string(split[0])) as temp_gff:
            fasta_input = SeqIO.to_dict(SeqIO.parse(temp_fasta, "fasta", generic_dna))
            gff_iter = GFF.parse(temp_gff, fasta_input)
            SeqIO.write(_check_gff(_fix_ncbi_id(gff_iter)), out_file, "genbank")
    gff_file.close()
    return

def run_panx(input_files, out_dir, ncpus=1, verbose=False):

    # create directory for input files
    input_file_dir = out_dir + "panx_run"
    if not os.path.exists(input_file_dir):
        os.mkdir(input_file_dir)

    # collect input files
    for f in input_files:
        convert_to_gbk(f, input_file_dir)

    # run panX
    cmd = ("panX.py" +
        " -fn " + input_file_dir +
        " -sl panx_run " +
        " -t " + str(ncpus)
        )

    if verbose:
        print("running cmd: ", cmd)

    subprocess.run(cmd, shell=True, check=True)

    return



def main():
    parser = argparse.ArgumentParser(description="""Runs PanX on GFF3 files and reformats output matrix.""")
    parser.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=str)

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="input GFF3 files (usually output from running Prokka)",
        type=str,
        nargs='+')

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    args = parser.parse_args()

    run_panx(args.input_files, args.output_dir, 
        ncpus=args.n_cpu, verbose=False)

if __name__ == '__main__':
    main()