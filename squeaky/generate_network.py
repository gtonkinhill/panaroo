from Bio import SeqIO
from collections import Counter, defaultdict
import networkx as nx
from merge_nodes import merge_nodes


def generate_network(cluster_file,
                     data_file,
                     prot_seq_file,
                     split_paralogs=True):

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
            if line[0] in seq_to_cluster:
                # this is a cluster centroid so keep it
                cluster_centroid_data[seq_to_cluster[line[0]]] = {
                    'prot_sequence': line[2],
                    'dna_sequence': line[3],
                    'annotation': line[4],
                    'description': line[5],
                }

    # load headers which contain adjacency information
    seq_ids = []
    with open(prot_seq_file, 'rU') as infile:
        for rec in SeqIO.parse(prot_seq_file, "fasta"):
            seq_ids.append(str(rec.id))

    # build graph using adjacency information and optionally split paralogs
    G = nx.Graph()
    n_nodes = len(cluster_members)
    temp_nodes = []
    for i, id in enumerate(seq_ids):
        current_cluster = seq_to_cluster[id]
        genome_id = id.split("_")[0]
        if id.split("_")[-1] == "0":
            # we're at the start of a contig
            prev = current_cluster
            if G.has_node(prev):
                G.node[prev]['size'] += 1
                G.node[prev]['members'].append(genome_id)
            else:
                if (prev in paralogs) and split_paralogs:
                    # create a new paralog
                    n_nodes += 1
                    prev = n_nodes
                    temp_nodes.append(prev)
                # add non paralog node
                G.add_node(
                    prev,
                    size=1,
                    centroid=cluster_centroids[current_cluster],
                    members=[genome_id],
                    protein=cluster_centroid_data[current_cluster]
                    ['prot_sequence'],
                    dna=cluster_centroid_data[current_cluster]['dna_sequence'],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    paralog=(current_cluster in paralogs))
        else:
            is_paralog = current_cluster in paralogs
            if is_paralog and split_paralogs:
                # create a new paralog
                n_nodes += 1
                neighbour = n_nodes
                temp_nodes.append(neighbour)
                G.add_node(
                    neighbour,
                    size=1,
                    centroid=cluster_centroids[current_cluster],
                    members=[genome_id],
                    protein=cluster_centroid_data[current_cluster]
                    ['prot_sequence'],
                    dna=cluster_centroid_data[current_cluster]['dna_sequence'],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    paralog=True)
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
                        protein=cluster_centroid_data[current_cluster]
                        ['prot_sequence'],
                        dna=cluster_centroid_data[current_cluster]
                        ['dna_sequence'],
                        annotation=cluster_centroid_data[current_cluster]
                        ['annotation'],
                        description=cluster_centroid_data[current_cluster]
                        ['description'],
                        paralog=is_paralog)
                    # add edge between nodes
                    G.add_edge(
                        prev, current_cluster, weight=1, members=[genome_id])
                else:
                    G.node[current_cluster]['size'] += 1
                    G.node[current_cluster]['members'].append(genome_id)
                    if G.has_edge(prev, current_cluster):
                        G[prev][current_cluster]['weight'] += 1
                        G[prev][current_cluster]['members'].append(genome_id)
                    else:
                        G.add_edge(
                            prev,
                            current_cluster,
                            weight=1,
                            members=[genome_id])
                prev = current_cluster

    processed = set()

    if split_paralogs:
        # iterate through nodes and merge paralogs
        for node in temp_nodes:
            if node not in processed:
                # Merge neighbouring nodes if their centroids and edges match
                paths = nx.single_source_shortest_path_length(
                    G, source=node, cutoff=2)
                merge_list = []
                n_neigh = set(
                    [G.node[n]['centroid'] for n in G.neighbors(node)])
                for (v, l) in paths.items():
                    if (l == 2) and (
                            G.node[v]['centroid'] == G.node[node]['centroid']):
                        v_neigh = set(
                            [G.node[n]['centroid'] for n in G.neighbors(v)])
                        if v_neigh == n_neigh:  # edges match by centroid
                            merge_list.append(v)
                processed.add(node)
                for m in merge_list:
                    n_nodes += 1
                    processed.add(m)
                    merge_nodes(G, m, node, n_nodes)
                    node = n_nodes

    return G
