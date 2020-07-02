import os
import tempfile
import argparse

import networkx as nx
from tqdm import tqdm
from collections import defaultdict, Counter
import itertools

from panaroo.cdhit import run_cdhit_est

def generate_db(gene_data_file, graph_file, outdir, min_support, ncpu=1, quiet=False):

    G = nx.read_gml(graph_file)

    # run cd-hit-est on centroids to ensure we collapse paralogs with length variation
    temp_dir = outdir + "tmpxjc585m7/"#os.path.join(tempfile.mkdtemp(dir=outdir), "")
    # all_centroids = set()
    # for n in G.nodes():
    #     all_centroids |= set(G.nodes[n]['centroid'].split(";"))
    # with open(temp_dir + "centroids.fasta", 'w') as outfile, \
    #     open(gene_data_file, 'r') as infile:
    #     for line in infile:
    #         line = line.strip().split(',')
    #         if line[2] in all_centroids:
    #             outfile.write('>' + line[2] + '\n' +
    #                 line[5] + '\n')

    # run_cdhit_est(
    #     input_file=temp_dir + "centroids.fasta",
    #     output_file=temp_dir + "cdhit_centroids",
    #     id=0.95,
    #     n_cpu=ncpu,
    #     quiet=quiet)

    # process the output
    clusters = []
    with open(temp_dir + "cdhit_centroids.clstr", 'r') as infile:
        c = []
        for line in infile:
            if line[0] == ">":
                clusters.append(c)
                c = []
            else:
                c.append(line.split(">")[1].split("...")[0])
        clusters.append(c)
    clusters = clusters[1:]

    # iterate through genes and group centroids that appear together
    tempG = nx.Graph()
    # add edges from cdhit clusters
    for cluster in clusters:
        for nA, nB in itertools.combinations(cluster, 2):
                tempG.add_edge(nA, nB)

    centroids_to_genes = defaultdict(set)
    centroids_to_description = defaultdict(set)
    for n in G.nodes():
        if int(G.nodes[n]['size']) < min_support: continue
        centroids = G.nodes[n]['centroid'].split(";")
        if len(centroids) <= 1:
            tempG.add_node(centroids[0])
        else:
            for nA, nB in itertools.combinations(centroids, 2):
                tempG.add_edge(nA, nB)
        for sid in centroids:
            centroids_to_genes[sid].add(G.nodes[n]['name'])
            for d in G.nodes[n]['description'].split(';'):
                centroids_to_description[sid].add(d)
    clusters = [list(comp) for comp in nx.connected_components(tempG)]

    # name clustes based on genes
    centroids_to_gene_name = {}
    centroids_to_gene_id = {}
    for gid, cluster in enumerate(clusters):
        name = set()
        for sid in cluster:
            name |= centroids_to_genes[sid]
        name = ";".join(sorted(list(name)))
        for sid in cluster:
            centroids_to_gene_name[sid] = name
            centroids_to_gene_id[sid] = str(gid)

    # run through gene_data and pull out the sequences in ariba format
    # with open(outdir+ 'ariba_db.fa', 'w') as ariba_fasta, \
    #     open(outdir + 'ariba_meta.tsv', 'w') as ariba_meta, \
    #     open(gene_data_file, 'r') as infile:
    #     for line in infile:
    #         line = line.strip().split(',')
    #         if line[2] in centroids_to_gene_name:
    #             seqname = line[3] + ';' + line[2]
    #             ariba_fasta.write('>' + seqname + '\n')
    #             ariba_fasta.write(line[5] + '\n')

    #             ariba_meta.write('\t'.join([seqname, '1', '0', '.', 
    #                 centroids_to_gene_name[line[2]], 
    #                 ';'.join(centroids_to_description[line[2]])]) + '\n')


    # run through gene_data and pull out the sequences in srst2
    # [clusterUniqueIdentifier]__[clusterSymbol]__[alleleSymbol]__[alleleUniqueIdentifier]
    gene_count = Counter()
    with open(outdir+ 'srst2_db.fa', 'w') as srst2_fasta, \
        open(gene_data_file, 'r') as infile:
        for line in infile:
            line = line.strip().split(',')
            if line[2] in centroids_to_gene_name:
                seqname = line[3] + ';' + line[2]
                srst2_fasta.write('>' + "__".join([
                    centroids_to_gene_id[line[2]], 
                    centroids_to_gene_name[line[2]],
                    seqname, 
                    str(gene_count[centroids_to_gene_id[line[2]]])]) + '\n')
                srst2_fasta.write(line[5] + '\n')
                gene_count[centroids_to_gene_id[line[2]]] += 1

    print("Final database includes", len(gene_count), 'families')

    return


def get_options():
    import argparse

    description = 'Create an Ariba database from a Panaroo pangenome'
    parser = argparse.ArgumentParser(description=description,
                                     prog='generate_ariba_db')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-d",
        "--directory",
        dest="directory",
        required=True,
        help="Location of Panaroo output directory",
        type=str)

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of a new output directory",
                         type=str)

    io_opts.add_argument("--min_support",
                         dest="min_support",
                         help="minimum number of genomes supporting a cluster for it to be included in the database (default=1)",
                         type=int,
                         default=1)

    # Other options
    parser.add_argument("-t",
                    "--threads",
                    dest="n_cpu",
                    help="number of threads to use (default=1)",
                    type=int,
                    default=1)

    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    args.directory = os.path.join(args.directory, "")

    # check files exist
    gene_data_file = args.directory + 'gene_data.csv'
    graph_file = args.directory + 'final_graph.gml'

    if not os.path.isfile(gene_data_file):
        raise RuntimeError("Missing gene_data.csv file!")
    if not os.path.isfile(graph_file):
        raise RuntimeError("Missing final_graph.gml file!")

    # generate databse
    generate_db(gene_data_file=gene_data_file,
                graph_file=graph_file,
                outdir=args.output_dir,
                min_support=args.min_support,
                ncpu=args.n_cpu,
                quiet=args.quiet)
    

    return


if __name__ == '__main__':
    main()
