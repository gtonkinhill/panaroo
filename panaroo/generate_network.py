from Bio import SeqIO
from collections import Counter, defaultdict
import networkx as nx
from panaroo.merge_nodes import merge_nodes
from panaroo.clean_network import collapse_paralogs
import numpy as np
from scipy.sparse import csc_matrix, lil_matrix

def generate_network(cluster_file, data_file, prot_seq_file, all_dna=False):

    # associate sequences with their clusters
    seq_to_cluster = {}
    cluster_centroids = {}
    cluster_members = defaultdict(list)
    with open(cluster_file, 'rU') as infile:
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
    with open(data_file, 'rU') as infile:
        next(infile)  # skip header
        for line in infile:
            line = line.strip().split(",")
            if line[2] in seq_to_cluster:
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
    current_context = None
    prev = None

    for i, id in enumerate(seq_ids):
        current_cluster = seq_to_cluster[id]
        loc = id.split("_")
        genome_id = loc[0]
        if loc[-1] == "0":

            if current_context is not None:
                for para, index in current_paralogs:
                    context_array = lil_matrix((1, n_centroids), dtype=np.int8)
                    for t in current_context:
                        context_array[0,centroid_index[t[0]]] = abs(index-t[1])
                    context_array = csc_matrix(context_array)
                    centroid_context[para[0]].append([para[1], para[2], 
                        para[3], context_array])
            current_paralogs = []
            current_context = [(cluster_centroids[current_cluster], int(loc[-1]))]
            
            # we're at the start of a contig
            if prev is not None: G.node[prev]['hasEnd'] = True
            prev = current_cluster
            if G.has_node(prev) and (prev not in paralogs):
                G.node[prev]['size'] += 1
                G.node[prev]['members'].append(genome_id)
                G.node[prev]['seqIDs'].append(id)
                G.node[prev]['hasEnd'] = True
                G.node[current_cluster]['lengths'].append(len(cluster_centroid_data[prev]
                        ['dna_sequence']))
                if all_dna:
                    G.node[prev]['dna'] += ";" + cluster_centroid_data[prev]['dna_sequence']
            else:
                if prev in paralogs:
                    # create a new paralog
                    n_nodes += 1
                    prev = n_nodes
                    temp_nodes.append(prev)
                    current_paralogs.append([[cluster_centroids[current_cluster], 
                        prev, genome_id, loc[1]], int(loc[2])])
                # add non paralog node
                G.add_node(
                    prev,
                    size=1,
                    centroid=cluster_centroids[current_cluster],
                    members=[genome_id],
                    seqIDs=[id],
                    hasEnd = True,
                    protein=cluster_centroid_data[current_cluster]
                    ['prot_sequence'],
                    dna=cluster_centroid_data[current_cluster]['dna_sequence'],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    lengths=[len(cluster_centroid_data[current_cluster]
                        ['dna_sequence'])],
                    paralog=(current_cluster in paralogs),
                    mergedDNA=False)
        else:
            current_context.append((cluster_centroids[current_cluster], int(loc[-1])))
            is_paralog = current_cluster in paralogs
            if is_paralog:
                # create a new paralog
                n_nodes += 1
                neighbour = n_nodes
                temp_nodes.append(neighbour)
                current_paralogs.append([[cluster_centroids[current_cluster], 
                    neighbour, genome_id, loc[1]], int(loc[2])])
                G.add_node(
                    neighbour,
                    size=1,
                    centroid=cluster_centroids[current_cluster],
                    members=[genome_id],
                    seqIDs=[id],
                    hasEnd = False,
                    protein=cluster_centroid_data[current_cluster]
                    ['prot_sequence'],
                    dna=cluster_centroid_data[current_cluster]
                    ['dna_sequence'],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    lengths=[len(cluster_centroid_data[current_cluster]
                        ['dna_sequence'])],
                    paralog=True,
                    mergedDNA=False)
                # add edge between nodes
                G.add_edge(prev, neighbour, weight=1, members=[genome_id])
                prev = neighbour
            else:
                if not G.has_node(current_cluster):
                    # we need to add the gene in
                    G.add_node(
                        current_cluster,
                        size=1,
                        centroid=cluster_centroids[current_cluster],
                        members=[genome_id],
                        seqIDs=[id],
                        hasEnd = False,
                        protein=cluster_centroid_data[current_cluster]
                        ['prot_sequence'],
                        dna=cluster_centroid_data[current_cluster]
                        ['dna_sequence'],
                        annotation=cluster_centroid_data[current_cluster]
                        ['annotation'],
                        description=cluster_centroid_data[current_cluster]
                        ['description'],
                        lengths=[len(cluster_centroid_data[current_cluster]
                        ['dna_sequence'])],
                        paralog=is_paralog,
                        mergedDNA=False)
                    # add edge between nodes
                    G.add_edge(prev,
                               current_cluster,
                               weight=1,
                               members=[genome_id])
                else:
                    G.node[current_cluster]['size'] += 1
                    G.node[current_cluster]['members'].append(genome_id)
                    G.node[current_cluster]['seqIDs'].append(id)
                    G.node[current_cluster]['lengths'].append(len(cluster_centroid_data[current_cluster]
                        ['dna_sequence']))
                    if all_dna:
                        G.node[current_cluster]['dna'] += ";" + cluster_centroid_data[current_cluster]['dna_sequence']
                    if G.has_edge(prev, current_cluster):
                        G[prev][current_cluster]['weight'] += 1
                        G[prev][current_cluster]['members'].append(genome_id)
                    else:
                        G.add_edge(prev,
                                   current_cluster,
                                   weight=1,
                                   members=[genome_id])
                prev = current_cluster

    return G, centroid_context
