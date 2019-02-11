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


    # build graph using adjacency information and optionally split paralogs
    G = nx.Graph()
    n_nodes = len(cluster_members)
    with open(prot_seq_file, 'rU') as infile:
        for rec in SeqIO.parse(prot_seq_file, "fasta"):
            id = str(rec.id)
            assembly_id = id.split("_")[0]

            if id.split("_")[-1]=="1":
                current_cluster = seq_to_cluster[id]
                prev = current_cluster
                if not G.has_node(prev):
                    if prev in paralogs:
                        prev = n_nodes + 1
                        n_nodes+=1
                    G.add_node(prev,
                        size=1,
                        centroid=cluster_centroids[current_cluster],
                        members=[id],
                        protein=cluster_centroid_prot_seq[current_cluster],
                        dna=cluster_centroid_dna_seq[current_cluster],
                        paralog=(current_cluster in paralogs))
                else:
                    G.node[prev]['size'] += 1
                    G.node[prev]['members'].append(id)
            else:
                current_cluster = seq_to_cluster[id]
                is_paralog = current_cluster in paralogs
                if is_paralog and split_paralogs:
                    # check if previous node neighbours include paralog
                    found = False
                    for neighbour in G.neighbors(prev):
                        if G.node[neighbour]["centroid"]==cluster_centroids[current_cluster]:
                            # we've found a previous instance of this paralog
                            found = True
                            break
                    if found:
                        G.node[neighbour]['size'] += 1
                        G.node[neighbour]['members'].append(id)
                        G[prev][neighbour]['members'].append(id)
                        G[prev][neighbour]['weight'] += 1
                    else:
                        # create new instance of the paralog
                        n_nodes+=1
                        neighbour = n_nodes
                        G.add_node(neighbour,
                            size=1,
                            centroid=cluster_centroids[current_cluster],
                            members=[id],
                            protein=cluster_centroid_prot_seq[current_cluster],
                            dna=cluster_centroid_dna_seq[current_cluster],
                            paralog=True)
                        # add edge between nodes
                        G.add_edge(prev, neighbour,
                            weight=1, members=[id])
                    prev = neighbour
                else:
                    if not G.has_node(current_cluster):
                        # we need to add the gene in
                        G.add_node(current_cluster,
                            size=1,
                            centroid=cluster_centroids[current_cluster],
                            members=[id],
                            protein=cluster_centroid_prot_seq[current_cluster],
                            dna=cluster_centroid_dna_seq[current_cluster],
                            paralog=is_paralog)
                        # add edge between nodes
                        G.add_edge(prev, current_cluster,
                            weight=1, members=[id])
                    else:
                        G.node[current_cluster]['size'] += 1
                        G.node[current_cluster]['members'].append(id)
                        if G.has_edge(prev, current_cluster):
                            G[prev][current_cluster]['weight'] += 1
                            G[prev][current_cluster]['members'].append(id)
                        else:
                            G.add_edge(prev, current_cluster,
                                weight=1, members=[id])
                    prev = current_cluster

    # fix additional trailing due to paralog
    if split_paralogs:
        bad_nodes = []
        for node in G.nodes():
            centroids = []
            members = []
            for neighbor in G.neighbors(node):
                centroids.append(G.node[neighbor]['centroid'])
                members.append(G.node[neighbor]['members'])
            if (len(set(centroids))<=1) and (len(centroids)>1):
                G.node[neighbor]['size'] = len(centroids)
                G.node[neighbor]['members'] = members
                G[node][neighbor]['weight'] = len(centroids)
                G[node][neighbor]['members'] = members
                for n in G.neighbors(node):
                    if n!=neighbor:
                        bad_nodes.append(n)
        for node in bad_nodes:
            G.remove_node(node)


    return G
