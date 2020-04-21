import os
import tempfile
import shutil
import argparse

import networkx as nx
from tqdm import tqdm
from joblib import Parallel, delayed
from collections import defaultdict, Counter
import itertools
import math
import numpy as np
from intbitset import intbitset

from .isvalid import *
from .__init__ import __version__
from .cdhit import run_cdhit
from .generate_output import *
from .clean_network import *
from .merge_nodes import merge_node_cluster, gen_edge_iterables, gen_node_iterables, iter_del_dups, del_dups


def make_list(inp):
    if not type(inp) == list:
        inp = [inp]
    return inp


def update_sid(sid, member_count):
    sid = sid.split("_")
    sid[0] = str(member_count + int(sid[0]))
    return ("_".join(sid))


def load_graphs(graph_files, n_cpu=1):
    for graph_file in graph_files:
        if not os.path.isfile(graph_file):
            print("Missing:", graph_file)
            raise RuntimeError("Missing graph file!")

    graphs = [nx.read_gml(graph_file) for graph_file in tqdm(graph_files)]
    isolate_names = list(
        itertools.chain.from_iterable(
            [G.graph['isolateNames'] for G in graphs]))

    member_count = 0
    node_count = 0
    id_mapping = []
    for i, G in enumerate(graphs):
        id_mapping.append({})
        # relabel nodes to be consecutive integers from 1
        mapping = {}
        for n in G.nodes():
            mapping[n] = node_count
            node_count += 1
        G = nx.relabel_nodes(G, mapping, copy=True)

        # set up edge members and remove conflicts.
        for e in G.edges():
            G[e[0]][e[1]]['members'] = intbitset([
                m + member_count for m in conv_list(G[e[0]][e[1]]['members'])
            ])

        # set up node parameters and remove conflicts.
        max_mem = -1
        for n in G.nodes():
            G.nodes[n]['centroid'] = [
                update_sid(sid, member_count)
                for sid in G.nodes[n]['centroid'].split(";")
            ]
            new_ids = set()
            for sid in conv_list(G.nodes[n]['seqIDs']):
                nid = update_sid(sid, member_count)
                id_mapping[i][sid] = nid
                new_ids.add(nid)
            G.nodes[n]['seqIDs'] = new_ids
            G.nodes[n]['protein'] = del_dups(G.nodes[n]['protein'].replace(
                '*', 'J').split(";"))
            G.nodes[n]['dna'] = del_dups(G.nodes[n]['dna'].split(";"))
            G.nodes[n]['lengths'] = conv_list(G.nodes[n]['lengths'])
            G.nodes[n]['longCentroidID'][1] = update_sid(
                G.nodes[n]['longCentroidID'][1], member_count)
            G.nodes[n]['members'] = intbitset(
                [m + member_count for m in conv_list(G.nodes[n]['members'])])
            max_mem = max(max_mem, max(G.nodes[n]['members']))

        member_count = max_mem + 1
        graphs[i] = G

    return graphs, isolate_names, id_mapping


def cluster_centroids(graphs,
                      outdir,
                      directories,
                      id_mapping,
                      len_dif_percent=0.95,
                      identity_threshold=0.98,
                      n_cpu=1):

    # create the files we will need
    temp_input_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_input_file.close()
    temp_output_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_output_file.close()

    # create input for cdhit
    orig_ids = {}
    ids_len_stop = {}
    with open(temp_input_file.name, 'w') as outfile:
        for i, d in enumerate(directories):
            with open(d + "gene_data.csv", 'r') as infile:
                next(infile)
                for line in infile:
                    line = line.split(",")
                    if line[2] not in id_mapping[i]:
                        continue  #its been filtered
                    orig_ids[id_mapping[i][line[2]]] = line[3]
                    ids_len_stop[id_mapping[i][line[2]]] = (len(
                        line[4]), "*" in line[4][1:-3])
                    outfile.write(">" + id_mapping[i][line[2]] + "\n" +
                                  line[4].replace("*", "J") + "\n")

    # Run cd-hit
    run_cdhit(temp_input_file.name,
              temp_output_file.name,
              id=identity_threshold,
              s=len_dif_percent,
              accurate=True,
              min_length=5,
              n_cpu=n_cpu)

    # Process output
    clusters = []
    with open(temp_output_file.name + ".clstr", 'r') as infile:
        c = []
        for line in infile:
            if line[0] == ">":
                clusters.append(c)
                c = []
            else:
                centroid = line.split(">")[1].split("...")[0]
                c.append(centroid)
        clusters.append(c)
    clusters = clusters[1:]

    # remove temporary files
    os.remove(temp_input_file.name)
    os.remove(temp_output_file.name)
    os.remove(temp_output_file.name + ".clstr")

    # rename centroids
    seqid_to_centroid = {}
    centroids_to_nodes = defaultdict(list)
    for cluster in clusters:
        for sid in cluster:
            seqid_to_centroid[sid] = cluster[0]
    all_centroids = set(list(seqid_to_centroid.values()))

    centroid_to_seqs = {}
    for i, d in enumerate(directories):
        with open(d + "gene_data.csv", 'r') as infile:
            next(infile)
            for line in infile:
                line = line.split(",")
                if line[2] not in id_mapping[i]: continue  #its been filtered
                if id_mapping[i][line[2]] in all_centroids:
                    centroid_to_seqs[id_mapping[i][line[2]]] = (line[4],
                                                                line[5])

    for G in graphs:
        for node in G.nodes():
            G.nodes[node]["centroid"] = list(
                set([
                    seqid_to_centroid[sid] for sid in G.nodes[node]['seqIDs']
                ]))
            G.nodes[node]["dna"] = [
                centroid_to_seqs[sid][1] for sid in G.nodes[node]["centroid"]
            ]
            G.nodes[node]["protein"] = [
                centroid_to_seqs[sid][0] for sid in G.nodes[node]["centroid"]
            ]
            G.nodes[node]["longCentroidID"] = max([
                (len(seq), sid) for seq, sid in zip(G.nodes[node]["dna"],
                                                    G.nodes[node]["centroid"])
            ])

    # determine node clusters
    tempG = nx.Graph()
    for centroid in centroids_to_nodes:
        if len(centroids_to_nodes[centroid]) <= 1:
            tempG.add_node(centroids_to_nodes[centroid][0])
        else:
            for nA, nB in itertools.combinations(centroids_to_nodes[centroid],
                                                 2):
                tempG.add_edge(nA, nB)

    clusters = [list(comp) for comp in nx.connected_components(tempG)]

    return clusters, seqid_to_centroid, orig_ids, ids_len_stop


def simple_merge_graphs(graphs, clusters):
    # Here, we only merge nodes that don't conflict

    # merge graphs
    merged_G = nx.compose_all(graphs)
    node_count = max(merged_G.nodes()) + 1

    for cluster in clusters:

        # check if there are any to collapse
        if len(cluster) <= 1: continue

        # check for conflicts
        seen = merged_G.nodes[cluster[0]]['members'].copy()
        noconflict = True
        for n in cluster[1:]:
            if not seen.isdisjoint(merged_G.nodes[n]['members']):
                noconflict = False
                break
            seen |= merged_G.nodes[n]['members']

        if noconflict:
            # no conflicts so merge
            node_count += 1
            merged_G = merge_node_cluster(merged_G,
                                          cluster,
                                          node_count,
                                          multi_centroid=True,
                                          check_merge_mems=True)

    return merged_G


def get_options():
    import argparse

    description = 'Merge independent runs of Panaroo'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo_merge_graphs')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-d",
        "--directories",
        dest="directories",
        required=True,
        help="Location of seperate Panaroo output directories",
        nargs='+')

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of a new output directory",
                         type=lambda x: is_valid_folder(parser, x))

    matching = parser.add_argument_group('Matching')

    matching.add_argument("-c",
                          "--threshold",
                          dest="id",
                          help="sequence identity threshold (default=0.95)",
                          default=0.98,
                          type=float)

    matching.add_argument(
        "-f",
        "--family_threshold",
        dest="family_threshold",
        help="protein family sequence identity threshold (default=0.7)",
        default=0.7,
        type=float)

    matching.add_argument("--len_dif_percent",
                          dest="len_dif_percent",
                          help="length difference cutoff (default=0.95)",
                          default=0.95,
                          type=float)

    matching.add_argument("--merge_paralogs",
                          dest="merge_paralogs",
                          help="don't split paralogs",
                          action='store_true',
                          default=False)

    matching.add_argument(
        "--length_outlier_support_proportion",
        dest="length_outlier_support_proportion",
        help=
        ("proportion of genomes supporting a gene with a length more " +
         "than 1.5x outside the interquatile range for genes in the same cluster"
         +
         " (default=0.01). Genes failing this test will be re-annotated at the "
         + "shorter length"),
        type=float,
        default=0.1)

    parser.add_argument(
        "--min_edge_support_sv",
        dest="min_edge_support_sv",
        help=(
            "minimum edge support required to call structural variants" +
            " in the presence/absence sv file (default=max(2, 0.01*n_samples))"
        ),
        default=2,
        type=int)

    # MSA options
    core = parser.add_argument_group('Gene alignment')
    core.add_argument(
        "-a",
        "--alignment",
        dest="aln",
        help=("Output alignments of core genes or all genes. Options are" +
              " 'core' and 'pan'. Default: 'None'"),
        type=str,
        choices=['core', 'pan'],
        default=None)
    core.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'prank', 'clustal', and default: 'mafft'",
        type=str,
        choices=['prank', 'clustal', 'mafft'],
        default="mafft")
    core.add_argument("--core_threshold",
                      dest="core",
                      help="Core-genome sample threshold (default=0.95)",
                      type=float,
                      default=0.95)

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
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    args.directories = [os.path.join(d, "") for d in args.directories]

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    print(
        "Merging graphs is still under active development and may change frequently!"
    )

    # Load graphs
    print("Loading graphs...")
    graphs, isolate_names, id_mapping = load_graphs(
        [d + "final_graph.gml" for d in args.directories], n_cpu=args.n_cpu)

    # cluster centroids
    print("Clustering centroids...")
    clusters, seqid_to_centroid, orig_ids, ids_len_stop = cluster_centroids(
        graphs=graphs,
        outdir=temp_dir,
        directories=args.directories,
        id_mapping=id_mapping,
        len_dif_percent=args.len_dif_percent,
        identity_threshold=args.id,
        n_cpu=args.n_cpu)

    # perform initial merge
    print("Performing inital merge...")
    G = simple_merge_graphs(graphs, clusters)

    print("Number of nodes in merged graph: ", G.number_of_nodes())

    print("Collapsing at DNA...")
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          dna_error_threshold=0.98,
                          correct_mistranslations=True,
                          length_outlier_support_proportion=args.
                          length_outlier_support_proportion,
                          n_cpu=args.n_cpu,
                          quiet=args.quiet,
                          depths = [1, 2])[0]

    print("Collapsing at families...")
    G = collapse_families(G,
                          seqid_to_centroid=seqid_to_centroid,
                          outdir=temp_dir,
                          family_threshold=args.family_threshold,
                          correct_mistranslations=False,
                          length_outlier_support_proportion=args.
                          length_outlier_support_proportion,
                          n_cpu=args.n_cpu,
                          quiet=args.quiet,
                          depths = [1, 2])[0]

    print("Number of nodes in merged graph: ", G.number_of_nodes())

    # Generate output
    print("Generating output...")

    # if requested merge paralogs
    if args.merge_paralogs:
        G = merge_paralogs(G)

    G.graph['isolateNames'] = isolate_names
    mems_to_isolates = {}
    for i, iso in enumerate(isolate_names):
        mems_to_isolates[i] = iso

    if not args.quiet:
        print("writing output...")

    # write out roary like gene_presence_absence.csv
    # get original annotaiton IDs, lengths and whether or
    # not an internal stop codon is present
    orig_ids = {}
    ids_len_stop = {}
    for i, d in enumerate(args.directories):
        with open(d + "gene_data.csv", 'r') as infile:
            next(infile)
            for line in infile:
                line = line.split(",")
                if line[2] not in id_mapping[i]: continue  #its been filtered
                orig_ids[id_mapping[i][line[2]]] = line[3]
                ids_len_stop[id_mapping[i][line[2]]] = (len(line[4]),
                                                        "*" in line[4][1:-3])

    G = generate_roary_gene_presence_absence(G,
                                             mems_to_isolates=mems_to_isolates,
                                             orig_ids=orig_ids,
                                             ids_len_stop=ids_len_stop,
                                             output_dir=args.output_dir)
    #Write out presence_absence summary
    generate_summary_stats(output_dir=args.output_dir)

    # write pan genome reference fasta file
    generate_pan_genome_reference(G,
                                  output_dir=args.output_dir,
                                  split_paralogs=False)

    # write out common structural differences in a matrix format
    generate_common_struct_presence_absence(
        G,
        output_dir=args.output_dir,
        mems_to_isolates=mems_to_isolates,
        min_variant_support=args.min_edge_support_sv)

    # add helpful attributes and write out graph in GML format
    for node in G.nodes():
        G.nodes[node]['size'] = len(G.nodes[node]['members'])
        G.nodes[node]['centroid'] = ";".join(G.nodes[node]['centroid'])
        G.nodes[node]['dna'] = ";".join(conv_list(G.nodes[node]['dna']))
        G.nodes[node]['protein'] = ";".join(conv_list(
            G.nodes[node]['protein']))
        G.nodes[node]['genomeIDs'] = ";".join(
            [str(m) for m in G.nodes[node]['members']])
        G.nodes[node]['geneIDs'] = ";".join(G.nodes[node]['seqIDs'])
        G.nodes[node]['degrees'] = G.degree[node]
        G.nodes[node]['members'] = list(G.nodes[node]['members'])
        G.nodes[node]['seqIDs'] = list(G.nodes[node]['seqIDs'])

    for edge in G.edges():
        G.edges[edge[0], edge[1]]['genomeIDs'] = ";".join(
            [str(m) for m in G.edges[edge[0], edge[1]]['members']])
        G.edges[edge[0],
                edge[1]]['members'] = list(G.edges[edge[0],
                                                   edge[1]]['members'])

    nx.write_gml(G, args.output_dir + "final_graph.gml")

    # write out merged gene_data and combined_DNA_CDS files
    with open(args.output_dir + "gene_data.csv", 'w') as outdata, \
    open(args.output_dir + "combined_DNA_CDS.fasta", 'w') as outdna:
        for i, d in enumerate(args.directories):
            with open(d + "gene_data.csv", 'r') as infile:
                header = next(infile)
                if i == 0: outdata.write(header)
                for line in infile:
                    line = line.strip().split(",")
                    if line[2] not in id_mapping[i]:
                        continue  #its been filtered
                    line[2] = id_mapping[i][line[2]]
                    outdata.write(",".join(line) + "\n")
                    outdna.write(">" + line[2] + "\n")
                    outdna.write(line[5] + "\n")

    # #Write out core/pan-genome alignments
    if args.aln == "pan":
        if not args.quiet: print("generating pan genome MSAs...")
        generate_pan_genome_alignment(G, temp_dir, args.output_dir, args.n_cpu,
                                      args.alr, isolate_names)
        core_nodes = get_core_gene_nodes(G, args.core, len(isolate_names))
        concatenate_core_genome_alignments(core_nodes, args.output_dir)
    elif args.aln == "core":
        if not args.quiet: print("generating core genome MSAs...")
        generate_core_genome_alignment(G, temp_dir, args.output_dir,
                                       args.n_cpu, args.alr, isolate_names,
                                       args.core, len(isolate_names))

    return


if __name__ == '__main__':
    main()
