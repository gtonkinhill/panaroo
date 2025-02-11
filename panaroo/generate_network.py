from Bio import SeqIO
from collections import Counter, defaultdict
import networkx as nx
from panaroo.clean_network import collapse_paralogs
import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from intbitset import intbitset


def generate_network(cluster_file, data_file, prot_seq_file, all_dna=False):

    # associate sequences with their clusters
    seq_to_cluster = {}
    seqid_to_centroid = {}
    cluster_centroids = {}
    cluster_members = defaultdict(list)
    with open(cluster_file, 'r') as infile:
        for line in infile:
            if line[0] == ">":
                cluster = int(line.strip().split()[-1])
            else:
                seq = line.split(">")[1].split("...")[0]
                seq_to_cluster[seq] = cluster
                cluster_members[cluster].append(seq.split("_")[0])
                if line.strip().split()[-1] == "*":
                    cluster_centroids[cluster] = seq

    # determine paralogs if required
    paralogs = set()
    for clust in cluster_members:
        genomes = [s.split("_")[0] for s in cluster_members[clust]]
        if len(genomes) != len(set(genomes)):
            paralogs.add(clust)

    # Load meta data such as sequence and annotation
    cluster_centroid_data = {}
    centroid_ids = set(cluster_centroids.values())
    with open(data_file, 'r') as infile:
        next(infile)  # skip header
        for line in infile:
            line = line.strip().split(",")
            if line[2] in centroid_ids:
                # this is a cluster centroid so keep it
                cluster_centroid_data[seq_to_cluster[line[2]]] = {
                    'prot_sequence': line[4],
                    'dna_sequence': line[5],
                    'annotation': line[6],
                    'description': line[7],
                }

    # create a dictionary of indexes for paralog context
    centroid_index = {}
    for i, centroid in enumerate(cluster_centroids.values()):
        centroid_index[centroid] = i
    n_centroids = len(centroid_index.keys())

    # load headers which contain adjacency information
    seq_ids = []
    for rec in SeqIO.parse(prot_seq_file, "fasta"):
        seq_ids.append(str(rec.id))

    # build graph using adjacency information and optionally split paralogs
    G = nx.Graph()
    centroid_context = defaultdict(list)
    n_nodes = len(cluster_members)
    temp_nodes = []
    prev = None

    for i, id in enumerate(seq_ids):
        current_cluster = seq_to_cluster[id]
        seqid_to_centroid[id] = cluster_centroids[current_cluster]
        loc = id.split("_")
        genome_id = int(loc[0])
        if loc[-1] == "0":
            # we're at the start of a contig
            if prev is not None: G.nodes[prev]['hasEnd'] = True
            prev = current_cluster
            if G.has_node(prev) and (prev not in paralogs):
                G.nodes[prev]['size'] += 1
                G.nodes[prev]['members'].add(genome_id)
                G.nodes[prev]['seqIDs'].add(id)
                G.nodes[prev]['hasEnd'] = True
                G.nodes[current_cluster]['lengths'].append(
                    len(cluster_centroid_data[prev]['dna_sequence']))
                if all_dna:
                    G.nodes[prev]['dna'] += [
                        cluster_centroid_data[prev]['dna_sequence']
                    ]
            else:
                if prev in paralogs:
                    # create a new paralog
                    n_nodes += 1
                    prev = n_nodes
                    temp_nodes.append(prev)
                    centroid_context[
                        cluster_centroids[current_cluster]].append(
                            [prev, genome_id])
                # add non paralog node
                G.add_node(
                    prev,
                    size=1,
                    centroid=[cluster_centroids[current_cluster]],
                    maxLenId=0,
                    members=intbitset([genome_id]),
                    seqIDs=set([id]),
                    hasEnd=True,
                    protein=[
                        cluster_centroid_data[current_cluster]['prot_sequence']
                    ],
                    dna=[
                        cluster_centroid_data[current_cluster]['dna_sequence']
                    ],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    lengths=[
                        len(cluster_centroid_data[current_cluster]
                            ['dna_sequence'])
                    ],
                    longCentroidID=(len(cluster_centroid_data[current_cluster]
                                        ['dna_sequence']),
                                    cluster_centroids[current_cluster]),
                    paralog=(current_cluster in paralogs),
                    mergedDNA=False)
        else:
            is_paralog = current_cluster in paralogs
            if is_paralog:
                # create a new paralog
                n_nodes += 1
                neighbour = n_nodes
                temp_nodes.append(neighbour)
                centroid_context[cluster_centroids[current_cluster]].append(
                    [neighbour, genome_id])
                G.add_node(
                    neighbour,
                    size=1,
                    centroid=[cluster_centroids[current_cluster]],
                    maxLenId=0,
                    members=intbitset([genome_id]),
                    seqIDs=set([id]),
                    hasEnd=False,
                    protein=[
                        cluster_centroid_data[current_cluster]['prot_sequence']
                    ],
                    dna=[
                        cluster_centroid_data[current_cluster]['dna_sequence']
                    ],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    lengths=[
                        len(cluster_centroid_data[current_cluster]
                            ['dna_sequence'])
                    ],
                    longCentroidID=(len(cluster_centroid_data[current_cluster]
                                        ['dna_sequence']),
                                    cluster_centroids[current_cluster]),
                    paralog=True,
                    mergedDNA=False)
                # add edge between nodes
                G.add_edge(prev,
                           neighbour,
                           size=1,
                           members=intbitset([genome_id]))
                prev = neighbour
            else:
                if not G.has_node(current_cluster):
                    # we need to add the gene in
                    G.add_node(
                        current_cluster,
                        size=1,
                        centroid=[cluster_centroids[current_cluster]],
                        maxLenId=0,
                        members=intbitset([genome_id]),
                        seqIDs=set([id]),
                        hasEnd=False,
                        protein=[
                            cluster_centroid_data[current_cluster]
                            ['prot_sequence']
                        ],
                        dna=[
                            cluster_centroid_data[current_cluster]
                            ['dna_sequence']
                        ],
                        annotation=cluster_centroid_data[current_cluster]
                        ['annotation'],
                        description=cluster_centroid_data[current_cluster]
                        ['description'],
                        lengths=[
                            len(cluster_centroid_data[current_cluster]
                                ['dna_sequence'])
                        ],
                        longCentroidID=(len(
                            cluster_centroid_data[current_cluster]
                            ['dna_sequence']),
                                        cluster_centroids[current_cluster]),
                        paralog=is_paralog,
                        mergedDNA=False)
                    # add edge between nodes
                    G.add_edge(prev,
                               current_cluster,
                               size=1,
                               members=intbitset([genome_id]))
                else:
                    G.nodes[current_cluster]['size'] += 1
                    G.nodes[current_cluster]['members'].add(genome_id)
                    G.nodes[current_cluster]['seqIDs'].add(id)
                    G.nodes[current_cluster]['lengths'].append(
                        len(cluster_centroid_data[current_cluster]
                            ['dna_sequence']))
                    if all_dna:
                        G.nodes[current_cluster]['dna'] += [
                            cluster_centroid_data[current_cluster]
                            ['dna_sequence']
                        ]
                    if G.has_edge(prev, current_cluster):
                        G[prev][current_cluster]['size'] += 1
                        G[prev][current_cluster]['members'].add(genome_id)
                    else:
                        G.add_edge(prev,
                                   current_cluster,
                                   size=1,
                                   members=intbitset([genome_id]))
                prev = current_cluster

    if prev is not None: 
        G.nodes[prev]['hasEnd'] = True

    return G, centroid_context, seqid_to_centroid
