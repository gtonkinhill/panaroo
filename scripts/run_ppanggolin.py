import argparse
import subprocess
import sys, os
import gffutils as gff

def run_ppanggolin(input_files, out_dir, defrag=False, ncpus=1, verbose=False):

    # create input file
    input_gff_locations = out_dir + "gff_file_locations.tab"
    with open(input_gff_locations, 'w') as outfile:
        for gff in input_files:
            prefix = os.path.splitext(os.path.basename(gff))[0]
            outfile.write(prefix + "\t" + gff + "\n")

    # run ppanggolin
    cmd = "ppanggolin workflow"
    cmd += " --anno " + input_gff_locations
    cmd += " -o " + out_dir
    cmd += " -c " + str(ncpus)
    cmd += " --force"

    if defrag:
        cmd += " --defrag"

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


def post_process_fmt(input_files, out_dir):

    # get mapping between GFF ids and ppanggolin ids
    id_mapping = {}
    for f in input_files:
        prefix = os.path.splitext(os.path.basename(f))[0]

        #Split file and parse
        gff_file = open(f, 'r')
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

        gene_count = 0
        for entry in parsed_gff.all_features(featuretype=()):
            if "CDS" not in entry.featuretype: continue
            ppanggolin_id = prefix + "_CDS_" + str(gene_count).zfill(4)
            id_mapping[ppanggolin_id] = entry.id
            gene_count += 1

        gff_file.close()
    
    id_mapping[''] = ''

    with open(out_dir + "ppanggolin_gene_presence_absence.csv", 'w') as outfile:
        with open(out_dir + "matrix.csv", 'r') as infile:
            outfile.write(next(infile))
            for line in infile:
                line = line.strip().split(",")
                for i in range(14, len(line)):
                    line[i] = ";".join([id_mapping[g.strip('"')] for g in line[i].split()])
                outfile.write(",".join(line) + "\n")

    return

def main():
    parser = argparse.ArgumentParser(
        description=
        """Runs ppanggolin on GFF3 files and reformats output matrix.""")
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

    parser.add_argument("--defrag",
                        dest="defrag",
                        help=("Turn ppanggolin defragmentation."),
                        action='store_true',
                        default=False)

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    args = parser.parse_args()

    args.output_dir = os.path.join(args.output_dir, "")

    run_ppanggolin(args.input_files,
                   args.output_dir,
                   defrag=args.defrag,
                   ncpus=args.n_cpu,
                   verbose=False)

    post_process_fmt(args.input_files, 
                     args.output_dir)

    return


if __name__ == '__main__':
    main()