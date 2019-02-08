from multi_prodigal import run_prodigal_multi
from cdhit import run_cdhit
from generate_network import generate_network
from clean_network import *
import os
import argparse
import tempfile
from Bio import SeqIO
import shutil
import networkx as nx


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg

def is_valid_folder(parser, arg):
    if not os.path.isdir(arg):
        parser.error("The folder %s does not exist!" % arg)
    else:
        return arg



def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--threshold", dest="id",
                    help="sequence identity threshold (default=0.95)",
                    type=float, default=0.95)

    parser.add_argument("-i", "--input", dest="input_files", required=True,
                    help="input files",
                    type=argparse.FileType('rU'), nargs='+')

    parser.add_argument("--train_dir", dest="train_dir",
                    help=("location of a training directory (a previous run" +
                        " of squeaky)"),
                    type=lambda x: is_valid_folder(parser, x))

    parser.add_argument("-o", "--out_dir", dest="output_dir", required=True,
                    help="location of an output directory",
                    type=lambda x: is_valid_folder(parser, x))

    parser.add_argument("--min_trailing_support",
                    dest="min_trailing_support",
                    help=("minimum cluster size to keep a gene called at the "
                        + "end of a contig (default=2)"),
                    type=int, default=2)

    parser.add_argument("--json",
                    dest="write_json",
                    help="write final graph out in JSON format.",
                    action='store_true',
                    default=False)

    parser.add_argument("-t", "--threads", dest="n_cpu",
                    help="number of threads to use (default=1)",
                    type=int, default=1)

    args = parser.parse_args()

    # Get prefix
    # prefix = os.path.splitext(os.path.basename(args.input_file))[0]

    # Create temporary directory
    temp_dir = tempfile.mkdtemp(dir=args.output_dir)

    # Rename sequences for simplicity later on and keep a record of it
    ref_file = args.output_dir + "/" + "sequence_name_reference.csv"

    # read in existing reference file if it exists
    prev_files = set()
    reference_list = []
    if os.path.exists(ref_file):
        with open(ref_file, 'rU') as infile:
            for line in infile:
                prev_files.add(line.split(",")[0])
                reference_list.append(line.strip().split(","))

    file_count = len(prev_files) - 1
    temp_input_files = []
    for f in args.input_files:
        temp_file = tempfile.NamedTemporaryFile(delete = False, dir=temp_dir,
                                                mode='w+t')
        contig_count=0
        file_count+=1
        for rec in SeqIO.parse(f, "fasta"):
            new_id = "_".join([str(file_count), str(contig_count)])
            temp_file.writelines([">" + new_id + "\n", str(rec.seq) + "\n"])
            reference_list.append((os.path.splitext(os.path.basename(
                f.name))[0], str(rec.id), file_count, contig_count))
            contig_count+=1
        temp_input_files.append(temp_file.name)
        temp_file.close()

    # write out updated reference file
    with open(ref_file, 'w') as outfile:
        for ref in reference_list:
            outfile.write(",".join(map(str,ref)) + "\n")

    # Run prodigal on new sequences
    run_prodigal_multi(temp_input_files, args.output_dir, temp_dir,
        n_cpu=args.n_cpu)

    # Cluster protein sequences using cdhit
    cd_hit_out = args.output_dir + "/" + "combined_protein_cdhit_out.txt"
    run_cdhit(input_file=args.output_dir + "combined_protein_CDS.fasta",
        output_file=cd_hit_out,
        id=args.id,
        n_cpu=args.n_cpu)

    # generate network from clusters and adjacency information
    G = generate_network(cd_hit_out+".clstr",
            args.output_dir + "combined_protein_CDS.fasta",
            args.output_dir + "combined_DNA_CDS.fasta")

    # write out raw edge list prior to filtering
    with open(args.output_dir + "/" + "network_edges_prefilter.csv", 'w') as outfile:
        outfile.write("source,target,count\n")
        for node1, node2, data in G.edges(data=True):
            outfile.write(",".join(map(str, [node1, node2, data['weight']])) +
                "\n")

    # remove low support trailing ends
    G = trim_low_support_trailing_ends(G, min_support=args.min_trailing_support,
        max_recursive=2)


    # write out final edge list and node attributes
    with open(args.output_dir + "/" + "network_edges.csv", 'w') as outfile:
        outfile.write("source,target,count\n")
        for node1, node2, data in G.edges(data=True):
            outfile.write(",".join(map(str, [node1, node2, data['weight']])) +
                "\n")

    with open(args.output_dir + "/" + "gene_cluster_attributes.csv", 'w') as outfile:
        outfile.write("Id,Label,size,centroid,protein,DNA\n")
        for node, data in G.nodes(data=True):
            outfile.write(",".join(map(str, [node,
                                             data['centroid'],
                                             data['size'],
                                             data['centroid'],
                                             data['protein'],
                                             data['dna']])) + "\n")

    # write out graph in GEXF format
    nx.write_gexf(G, args.output_dir + "/" + "final_graph.gexf")

    # # optionally write out graph in JSON format
    # if args.write_json:


    # remove temp TemporaryDirectory
    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
