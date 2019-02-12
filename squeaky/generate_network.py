from Bio import SeqIO
from collections import Counter, defaultdict
import networkx as nx


def generate_network(cluster_file, prot_seq_file, dna_seq_file,
    split_paralogs=True):

    # associate sequences with their clusters
    seq_to_cluster = {}
    cluster_centroids = {}
    cluster_members = defaultdict(list)
    with open(cluster_file, 'rU') as infile:
        for line in infile:
            if line[0]==">":
                cluster = int(line.strip().split()[-1])
            else:
                seq = line.split(">")[1].split("...")[0]
                seq_to_cluster[seq] = cluster
                cluster_members[cluster].append(seq.split("_")[0])
                if line.strip().split()[-1]=="*":
                    cluster_centroids[cluster] = seq

    # determine paralogs if required
    paralogs = set()
    for  clust in cluster_members:
        genomes = [s.split("_")[0]  for s in cluster_members[clust]]
        if len(genomes) != len(set(genomes)):
            paralogs.add(clust)

    # identify DNA seq centroids
    cluster_centroid_dna_seq = {}
    with open(dna_seq_file, 'rU') as infile:
        for rec in SeqIO.parse(infile, "fasta"):
            id = str(rec.id)
            if cluster_centroids[seq_to_cluster[id]]==id:
                cluster_centroid_dna_seq[seq_to_cluster[id]] = str(rec.seq)

    # identify protein seq centroids
    cluster_centroid_prot_seq = {}
    with open(prot_seq_file, 'rU') as infile:
        for rec in SeqIO.parse(infile, "fasta"):
            id = str(rec.id)
            if cluster_centroids[seq_to_cluster[id]]==id:
                cluster_centroid_prot_seq[seq_to_cluster[id]] = str(rec.seq)

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
        if id.split("_")[-1]=="1":
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
                G.add_node(prev,
                    size=1,
                    centroid=cluster_centroids[current_cluster],
                    members=[genome_id],
                    protein=cluster_centroid_prot_seq[current_cluster],
                    dna=cluster_centroid_dna_seq[current_cluster],
                    paralog=(current_cluster in paralogs))
        else:
            is_paralog = current_cluster in paralogs
            if is_paralog and split_paralogs:
                # create a new paralog
                n_nodes+=1
                neighbour = n_nodes
                temp_nodes.append(neighbour)
                G.add_node(neighbour,
                    size=1,
                    centroid=cluster_centroids[current_cluster],
                    members=[genome_id],
                    protein=cluster_centroid_prot_seq[current_cluster],
                    dna=cluster_centroid_dna_seq[current_cluster],
                    paralog=True)
                # add edge between nodes
                G.add_edge(prev, neighbour,
                    weight=1, members=[genome_id])
                prev = neighbour
            else:
                if not G.has_node(current_cluster):
                    # we need to add the gene in
                    G.add_node(current_cluster,
                        size=1,
                        centroid=cluster_centroids[current_cluster],
                        members=[genome_id],
                        protein=cluster_centroid_prot_seq[current_cluster],
                        dna=cluster_centroid_dna_seq[current_cluster],
                        paralog=is_paralog)
                    # add edge between nodes
                    G.add_edge(prev, current_cluster,
                        weight=1, members=[genome_id])
                else:
                    G.node[current_cluster]['size'] += 1
                    G.node[current_cluster]['members'].append(genome_id)
                    if G.has_edge(prev, current_cluster):
                        G[prev][current_cluster]['weight'] += 1
                        G[prev][current_cluster]['members'].append(genome_id)
                    else:
                        G.add_edge(prev, current_cluster,
                            weight=1, members=[genome_id])
                prev = current_cluster

    processed=set()

    if split_paralogs:
        # iterate through nodes and merge paralogs
        for node in temp_nodes:
            if node not in processed:
                # Find neighbouring nodes and merge them if the centroids AND edges match
                paths = nx.single_source_shortest_path_length(G ,source=node, cutoff=2)
                merge_list = []
                n_neigh = set([G.node[n]['centroid'] for n in G.neighbors(node)])
                for (v, l) in paths.items():
                    if (l==2) and (G.node[v]['centroid']==G.node[node]['centroid']):
                        v_neigh = set([G.node[n]['centroid'] for n in G.neighbors(v)])
                        if v_neigh==n_neigh: # edges match by centroid
                            merge_list.append(v)
                processed.add(node)
                for m in merge_list:
                    n_nodes+=1
                    processed.add(m)
                    merge_nodes(G, m, node, n_nodes)
                    node = n_nodes

    return G


def merge_nodes(G, nodeA, nodeB, newNode, multi_centroid=False):
    # First create a new node and combine the attributes
    if multi_centroid:
        G.add_node(newNode,
            size=G.node[nodeA]['size'] + G.node[nodeB]['size'],
            centroid=set(G.node[nodeA]['centroid'], G.node[nodeB]['centroid']),
            members=G.node[nodeA]['members'] + G.node[nodeB]['members'],
            protein=set(G.node[nodeA]['protein'], G.node[nodeB]['protein']),
            dna=set(G.node[nodeA]['dna'], G.node[nodeB]['dna']),
            paralog=(G.node[nodeA]['paralog'] or G.node[nodeB]['paralog']))
    else:
        G.add_node(newNode,
            size=G.node[nodeA]['size'] + G.node[nodeB]['size'],
            centroid=G.node[nodeA]['centroid'],
            members=G.node[nodeA]['members'] + G.node[nodeB]['members'],
            protein=G.node[nodeA]['protein'],
            dna=G.node[nodeA]['dna'],
            paralog=(G.node[nodeA]['paralog'] or G.node[nodeB]['paralog']))

    # Now iterate through neighbours of each node and add them to the new node
    neigboursB = list(G.neighbors(nodeB))
    neigboursA = list(G.neighbors(nodeA))
    for neighbor in neigboursA:
        if neighbor in neigboursB:
            G.add_edge(newNode, neighbor,
                weight=G[nodeA][neighbor]['weight'] + G[nodeB][neighbor]['weight'],
                members=G[nodeA][neighbor]['members'] + G[nodeB][neighbor]['members'])
            neigboursB.remove(neighbor)
        else:
            G.add_edge(newNode, neighbor,
                weight=G[nodeA][neighbor]['weight'],
                members=G[nodeA][neighbor]['members'])

    for neighbor in neigboursB:
        G.add_edge(newNode, neighbor,
            weight=G[nodeB][neighbor]['weight'],
            members=G[nodeB][neighbor]['members'])

    # remove old nodes from Graph
    G.remove_nodes_from([nodeA, nodeB])

    return G
