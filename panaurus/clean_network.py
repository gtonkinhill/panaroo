import networkx as nx
from panaurus.cdhit import cluster_nodes_cdhit
from panaurus.merge_nodes import merge_nodes
from collections import defaultdict
from panaurus.cdhit import is_valid
from itertools import chain


# Genes at the end of contigs are more likely to be false positives thus
# we can remove those with low support
def trim_low_support_trailing_ends(G, min_support=3, max_recursive=2):

    # fix trailing
    for iter in range(max_recursive):
        bad_nodes = []
        removed = False
        for (node, val) in G.degree():
            if val <= 1:  # trailing node
                if G.node[node]['size'] < min_support:
                    bad_nodes.append(node)
        for node in bad_nodes:
            G.remove_node(node)
            removed = True

        # fix trailing due to paralog
        bad_nodes = []
        for node in G.nodes():
            centroids = []
            for neighbor in G.neighbors(node):
                centroids.append(G.node[neighbor]['centroid'])
            if len(set(centroids)) <= 1:
                if len(centroids) < min_support:
                    bad_nodes.append(node)
        for node in bad_nodes:
            G.remove_node(node)
            removed = True

        if not removed: break

    return G


def collapse_families(G,
                      outdir,
                      family_threshold=0.7,
                      dna_error_threshold=0.99,
                      correct_mistranslations=False,
                      n_cpu=1,
                      quiet=False):

    node_count = max(list(G.nodes())) + 10

    if correct_mistranslations:
        depths = [3]
    else:
        depths = range(1,10)
        
    for d in depths:
        search_space = set(G.nodes())

        while len(search_space) > 0:
            # look for nodes to merge
            temp_node_list = list(search_space)
            removed_nodes = set()
            for node in temp_node_list:
                if node in removed_nodes: continue

                if G.degree[node] <= 2:
                    search_space.remove(node)
                    removed_nodes.add(node)
                    continue

                # find neighbouring nodes and cluster their centroid with cdhit
                neighbours = [
                    v for u, v in nx.bfs_edges(G, source=node, depth_limit=d)
                ]

                if correct_mistranslations:
                    clusters = cluster_nodes_cdhit(
                        G,
                        neighbours,
                        outdir,
                        id=dna_error_threshold,
                        n_cpu=n_cpu,
                        quiet=True,
                        dna=True,
                        use_local=False,
                        # aS=0.9,
                        prevent_para=False,
                        accurate=False)
                else:
                    clusters = cluster_nodes_cdhit(G,
                                                    neighbours,
                                                    outdir,
                                                    id=family_threshold,
                                                    n_cpu=n_cpu,
                                                    quiet=True,
                                                    dna=False,
                                                    prevent_para=False)

                # for cluster in cluster_dict.values():
                for cluster in clusters:
                    # check if there are any to collapse
                    if len(cluster) <= 1: continue

                    # check for conflicts
                    members = list(
                        chain.from_iterable(
                            [G.node[n]['members'] for n in cluster]))

                    if (len(members) == len(set(members))):
                        # no conflicts so merge
                        node_count += 1
                        for neig in cluster:
                            removed_nodes.add(neig)
                            if neig in search_space: search_space.remove(neig)
                        temp_c = cluster.copy()
                        G = merge_nodes(G, temp_c.pop(), temp_c.pop(),
                                        node_count)
                        while (len(temp_c) > 0):
                            G = merge_nodes(G, node_count, temp_c.pop(),
                                            node_count + 1)
                            node_count += 1
                        search_space.add(node_count)
                    else:
                        if correct_mistranslations:
                            # merge if the centroids don't conflict and the nodes are adjacent in the conflicting genome
                            # this corresponds to a mistranslation where one gene has been split into two in a subset of genomes

                            # work out which nodes each genome has
                            member_to_nodes = defaultdict(list)
                            for n in cluster:
                                for mem in G.node[n]['members']:
                                    member_to_nodes[mem].append(n)

                            should_merge = True
                            for mem in member_to_nodes:
                                if len(member_to_nodes[mem]) <= 1: continue
                                temp_centroids = [
                                    G.node[n]['centroid']
                                    for n in member_to_nodes[mem]
                                ]
                                if len(temp_centroids) != len(
                                        set(temp_centroids)):
                                    # matching centroids so dont merge
                                    should_merge = False
                                sub_G = G.subgraph(member_to_nodes[mem])
                                if not nx.is_connected(sub_G):
                                    should_merge = False

                            if should_merge:
                                node_count += 1
                                for neig in cluster:
                                    removed_nodes.add(neig)
                                    if neig in search_space:
                                        search_space.remove(neig)
                                temp_c = cluster.copy()
                                G = merge_nodes(G,
                                                temp_c.pop(),
                                                temp_c.pop(),
                                                node_count,
                                                check_merge_mems=False)
                                while (len(temp_c) > 0):
                                    G = merge_nodes(G,
                                                    node_count,
                                                    temp_c.pop(),
                                                    node_count + 1,
                                                    check_merge_mems=False)
                                    node_count += 1
                                search_space.add(node_count)
                        else:
                            # there is a conflict in the merge, check if we can split based on neighbours
                            node_neighbours = defaultdict(list)
                            for n in cluster:
                                node_neighbours[tuple(
                                    sorted(list(G.neighbors(n))))].append(n)
                            for sub_cluster in node_neighbours.values():
                                if len(sub_cluster) <= 1: continue
                                sub_mems = list(
                                    chain.from_iterable([
                                        G.node[n]['members']
                                        for n in sub_cluster
                                    ]))
                                if len(sub_mems) != len(set(sub_mems)):
                                    continue
                                node_count += 1
                                for neig in sub_cluster:
                                    removed_nodes.add(neig)
                                    if neig in search_space:
                                        search_space.remove(neig)
                                temp_c = sub_cluster.copy()
                                G = merge_nodes(G, temp_c.pop(), temp_c.pop(),
                                                node_count)
                                while (len(temp_c) > 0):
                                    G = merge_nodes(G, node_count,
                                                    temp_c.pop(),
                                                    node_count + 1)
                                    node_count += 1
                                search_space.add(node_count)

                search_space.remove(node)

    return G


def collapse_paralogs(G, quiet=False):

    node_count = max(list(G.nodes())) + 10

    was_merged = True
    while was_merged:
        was_merged = False
        search_space = set(G.nodes())
        while len(search_space) > 0:
            # look for nodes to merge
            temp_node_list = list(search_space)
            removed_nodes = set()
            for node in temp_node_list:
                if node in removed_nodes: continue

                neigbour_centroids = defaultdict(list)
                # find neighbours centroids
                for neigh in [
                        v
                        for u, v in nx.bfs_edges(G, source=node, depth_limit=3)
                ]:
                    neigbour_centroids[G.node[neigh]['centroid']].append(neigh)

                for centroid in neigbour_centroids:
                    # check if there are any to collapse
                    if (len(neigbour_centroids[centroid]) > 1):
                        # check for conflicts
                        genomes = []
                        for n in neigbour_centroids[centroid]:
                            genomes += G.node[n]['members']
                        if len(genomes) == len(set(genomes)):
                            # no conflicts in merge
                            was_merged = True
                            node_count += 1
                            # merge neighbours with this centroid
                            for neigh in neigbour_centroids[centroid]:
                                removed_nodes.add(neigh)
                                if neigh in search_space:
                                    search_space.remove(neigh)
                            temp_c = neigbour_centroids[centroid].copy()
                            G = merge_nodes(G, temp_c.pop(), temp_c.pop(),
                                            node_count)
                            while (len(temp_c) > 0):
                                G = merge_nodes(G, node_count, temp_c.pop(),
                                                node_count + 1)
                                node_count += 1
                            search_space.add(node_count)
                        else:
                            # possible conflict (try splitting by neighbours)
                            neigh_neighbours = defaultdict(list)
                            for n in neigbour_centroids[centroid]:
                                neigh_neighbours[tuple(sorted(
                                    G.neighbors(n)))].append(n)

                            for nn in neigh_neighbours:
                                if len(neigh_neighbours[nn]) < 2: continue
                                # we've found subsets that share the same neighbours so merge
                                was_merged = True
                                node_count += 1
                                # merge neighbours with this centroid
                                for neigh in neigh_neighbours[nn]:
                                    removed_nodes.add(neigh)
                                    if neigh in search_space:
                                        search_space.remove(neigh)
                                temp_c = neigh_neighbours[nn].copy()
                                G = merge_nodes(G, temp_c.pop(), temp_c.pop(),
                                                node_count)
                                while (len(temp_c) > 0):
                                    G = merge_nodes(G, node_count,
                                                    temp_c.pop(),
                                                    node_count + 1)
                                    node_count += 1
                                search_space.add(node_count)

                search_space.remove(node)

    return G


def merge_paralogs(G):

    node_count = max(list(G.nodes())) + 10

    # group paralog nodes by centroid
    paralog_centroid_dict = defaultdict(list)
    for node in G.nodes():
        if G.node[node]['paralog']:
            paralog_centroid_dict[G.node[node]['centroid']].append(node)

    # merge paralog nodes that share the same centroid
    for centroid in paralog_centroid_dict:
        node_count += 1
        temp_c = paralog_centroid_dict[centroid]
        G = merge_nodes(G,
                        temp_c.pop(),
                        temp_c.pop(),
                        node_count,
                        check_merge_mems=False)
        while (len(temp_c) > 0):
            G = merge_nodes(G,
                            node_count,
                            temp_c.pop(),
                            node_count + 1,
                            check_merge_mems=False)
            node_count += 1

    return (G)


def clean_misassembly_edges(G, threshold):

    bad_edges = []
    for edge in G.edges():
        if float(G.edges[edge]['weight']) < (threshold * 
                min(int(G.node[edge[0]]['size']), int(G.node[edge[1]]['size']))):
            bad_edges.append(edge)
    
    for edge in bad_edges:
        G.remove_edge(edge[0], edge[1])
    
    return (G)