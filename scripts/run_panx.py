import argparse
import subprocess
import sys, os

def clean_gbk(gff_file_name, outdir):
    prefix = os.path.splitext(os.path.basename(gff_file_name))[0]
    prefix = prefix.replace("-", "")
    out_file = outdir+prefix+".gbk"

    lines = []
    source_lenghts = []
    with open(gff_file_name, 'r') as infile:
        for line in infile:
            if "     source" in line:
                source_lenghts.append(line.strip().split()[1].split("..")[1])
    
    source_counter = 0
    locus_count = 0
    prev_loc = False
    with open(out_file, 'w') as outfile:
        with open(gff_file_name, 'r') as infile:
            for line in infile:
                if line[:5]=="LOCUS":
                    outfile.write("LOCUS       " +
                        "contig" + str(locus_count) +
                        "           " + source_lenghts[source_counter] +
                        " bp    DNA     linear CON 19-AUG-2019\n")
                    prev_loc = True
                    locus_count += 1
                    source_counter += 1
                elif prev_loc:
                    if "DEFINITION" not in line:
                        prev_loc = False
                        continue
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)

    return prefix


def run_panx(input_files, out_dir, ncpus=1, verbose=False):

    # create directory for input files
    input_file_dir = out_dir + "panx_run/"
    if not os.path.exists(input_file_dir):
        os.mkdir(input_file_dir)

    # collect input files
    new_prefixes = []
    for f in input_files:
        new_prefixes.append(clean_gbk(f, input_file_dir))

    # run panX
    cmd = ("~/software/pan-genome-analysis/panX.py" +
        " -fn " + input_file_dir +
        " -sl panx_run " +
        " -t " + str(ncpus)
        )

    if verbose:
        print("running cmd: ", cmd)

    subprocess.run(cmd, shell=True, check=True)

    return new_prefixes

def post_process_fmt(new_prefixes, panx_dir, out_dir):

    # generate an index for files
    n_samples = len(new_prefixes)
    index_dict = {}
    for i, p in enumerate(new_prefixes):
        index_dict[p] = i

    with open(out_dir + "panx_gene_presence_absence.csv", 'w') as outfile:
        with open(panx_dir + "allclusters_final.tsv", 'r') as infile:
            outfile.write(",".join(new_prefixes) + "\n")
            for line in infile:
                pa = n_samples*[""]
                cluster = line.strip().split()
                for g in cluster:
                    g = g.split("|")
                    pa[index_dict[g[0]]] = g[1]
                outfile.write(",".join(pa) + "\n")

    return



def main():
    parser = argparse.ArgumentParser(description="""Runs PanX on GBK files and reformats output matrix.""")
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

    new_prefixes = run_panx(args.input_files, args.output_dir, 
        ncpus=args.n_cpu, verbose=False)

    post_process_fmt(new_prefixes, 
        args.output_dir + "panx_run/", args.output_dir)

    return

if __name__ == '__main__':
    main()