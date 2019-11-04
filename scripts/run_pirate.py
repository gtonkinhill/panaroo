import argparse
import subprocess
import sys, os
from shutil import copyfile
import gffutils as gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_pirate(input_files, out_dir, ncpus=1, verbose=False):

    # create directory for input files
    input_file_dir = out_dir + "pirate_run/"
    if not os.path.exists(input_file_dir):
        os.mkdir(input_file_dir)

    # collect input files
    for f in input_files:
        base = os.path.basename(f)
        base = base.replace("-", "_")
        copyfile(f, input_file_dir + base)

    # run panX
    cmd = ("PIRATE" +
        " -i " + input_file_dir +
        " -o " + input_file_dir +
        " -t " + str(ncpus)
        )

    if verbose:
        print("running cmd: ", cmd)

    subprocess.run(cmd, shell=True, check=True)

    return

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

def post_process_fmt(input_files, pirate_dir, out_dir):

    file_id_maps = {}
    for f in input_files:
        #Split file and parse
        gff_file=open(f, 'r')
        lines = gff_file.read()
        split = lines.split('##FASTA')

        if len(split) != 2:
            print("Problem reading GFF3 file: ", gff_file.name)
            raise RuntimeError("Error reading prokka input!")

        parsed_gff = gff.create_db(clean_gff_string(split[0]),
                                dbfn=":memory:",
                                force=True,
                                keep_order=True,
                                from_string=True)

        gene_count=1
        count_to_id = {}
        for entry in parsed_gff.all_features(featuretype=()):
            if "CDS" not in entry.featuretype: continue
            count_to_id[gene_count] = entry.id
            gene_count += 1

        file_id_maps[os.path.splitext(os.path.basename(f))[0].replace("-", "_")] = count_to_id

        gff_file.close()

    with open(pirate_dir + "PIRATE.gene_families.ordered.tsv", 'r') as infile:
        with open(out_dir + "pirate_gene_presence_absence.csv", 'w') as outfile:
            line = next(infile)
            line = line.strip().split("\t")
            outfile.write("\t".join([line[0]] + line[22:]) + "\n")
            
            for line in infile:
                line = line.strip().split("\t")
                for i in range(22, len(line)):
                    if line[i]=="": continue
                    # currently take the first as we do in panaroo for this format. May want to redo.
                    line[i] = line[i].split(";")[0]
                    sample = "_".join(line[i].split("_")[:-1])
                    gene_num = int(line[i].split("_")[-1])
                    line[i] = file_id_maps[sample][gene_num]
                outfile.write("\t".join([line[0]] + line[22:]) + "\n")

    return



def main():
    parser = argparse.ArgumentParser(description="""Runs PIRATE on GFF3 files and reformats output matrix.""")
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
        help="input GBK files (usually output from running Prokka)",
        type=str,
        nargs='+')

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    args = parser.parse_args()

    args.output_dir = os.path.join(args.output_dir, "")

    run_pirate(args.input_files, args.output_dir, 
        ncpus=args.n_cpu, verbose=False)

    post_process_fmt(args.input_files, 
        args.output_dir + "pirate_run/", args.output_dir)

    return

if __name__ == '__main__':
    main()