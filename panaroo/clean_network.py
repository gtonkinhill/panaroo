import networkx as nx
from panaroo.cdhit import cluster_nodes_cdhit, iterative_cdhit
from panaroo.merge_nodes import merge_nodes
from collections import defaultdict
from panaroo.cdhit import is_valid
from itertools import chain


# Genes at the end of contigs are more likely to be false positives thus
# we can remove those with low support
def trim_low_support_trailing_ends(G, min_support=3, max_recursive=2):

    # fix trailing
    for i in range(max_recursive):
        bad_nodes = []
        removed = False
        for (node, val) in G.degree():
            if val <= 1:  # trailing node
                if G.node[node]['size'] < min_support:
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
        depths = [1,2]
    else:
        depths = [1,2]  # range(1,10)

    # precluster for speed
    if correct_mistranslations:
        clusters = iterative_cdhit(
            G,
            G.nodes(),
            outdir,
            thresholds=[1,0.99],
            n_cpu=n_cpu,
            quiet=True,
            dna=True,
            # aL=0.6,
            use_local=False,
            accurate=False)
    else:
        clusters = iterative_cdhit(G,
            G.nodes(),
            outdir,
            thresholds=[0.99,0.95,0.90,0.8,0.7],
            n_cpu=n_cpu,
            quiet=True,
            dna=False)
    cluster_dict = {}
    for c, clust in enumerate(clusters):
        for n in clust:
            cluster_dict[n] = c

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

                temp_clusters = defaultdict(list)
                for neigh in neighbours:
                    temp_clusters[cluster_dict[neigh]].append(neigh)
                temp_clusters = temp_clusters.values()

                # temp_clusters = [[item for sublist in temp_clusters for item in sublist]]
                
                clusters = []
                for cluster in temp_clusters:
                    if len(cluster)<2:
                        clusters.append([cluster])
                        continue
                    if correct_mistranslations:
                        clusters += cluster_nodes_cdhit(
                            G,
                            cluster,
                            outdir,
                            id=dna_error_threshold,
                            n_cpu=n_cpu,
                            quiet=True,
                            dna=True,
                            use_local=False,
                            prevent_para=False,
                            accurate=False)
                    else:
                        clusters += cluster_nodes_cdhit(G,
                            cluster,
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
                                        node_count, 
                                        multi_centroid=(not correct_mistranslations))
                        cluster_dict[node_count] = cluster_dict[cluster[0]]
                        while (len(temp_c) > 0):
                            G = merge_nodes(G, node_count, temp_c.pop(),
                                            node_count + 1,
                                            multi_centroid=(not correct_mistranslations))
                            node_count += 1
                            cluster_dict[node_count] = cluster_dict[cluster[0]]
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
                                                multi_centroid=False,
                                                check_merge_mems=False)
                                cluster_dict[node_count] = cluster_dict[cluster[0]]
                                while (len(temp_c) > 0):
                                    G = merge_nodes(G,
                                                    node_count,
                                                    temp_c.pop(),
                                                    node_count + 1,
                                                    multi_centroid=False,
                                                    check_merge_mems=False)
                                    node_count += 1
                                    cluster_dict[node_count] = cluster_dict[cluster[0]]
                                search_space.add(node_count)
                        else:
                            # there is a conflict in the merge, check if we can split based on neighbours
                            was_merged = True
                            already_merged = set()
                            while was_merged:
                                was_merged = False
                                pos_merges = []
                                for nA in cluster:
                                    if nA in already_merged: continue
                                    best_inter = -1
                                    for nB in cluster:
                                        if nA == nB: continue
                                        if len(
                                                set(G.node[nA]
                                                    ['members']).intersection(
                                                        set(G.node[nB]
                                                            ['members']))) > 0:
                                            continue
                                        temp_inter = len(
                                            set(G.neighbors(nA)).intersection(
                                                set(G.neighbors(nB))))
                                        # if temp_inter==0: continue
                                        if temp_inter > best_inter:
                                            best_inter = temp_inter
                                            best_merge = nB
                                    if best_inter == -1:
                                        # none left to merge with this node
                                        already_merged.add(nA)
                                    else:
                                        pos_merges.append(
                                            (best_inter, nA, best_merge))
                                if len(pos_merges) > 0:
                                    was_merged = True
                                    best_merge = max(pos_merges)
                                    node_count += 1
                                    G = merge_nodes(G, best_merge[1],
                                                    best_merge[2], node_count)
                                    cluster_dict[node_count] = cluster_dict[cluster[0]]
                                    if best_merge[1] in search_space:
                                        search_space.remove(best_merge[1])
                                    if best_merge[2] in search_space:
                                        search_space.remove(best_merge[2])
                                    removed_nodes.add(best_merge[1])
                                    removed_nodes.add(best_merge[2])
                                    cluster.remove(best_merge[1])
                                    cluster.remove(best_merge[2])
                                    cluster.append(node_count)
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
                                neighs = tuple(sorted(G.neighbors(n)))
                                if len(neighs)>0:
                                    neigh_neighbours[neighs].append(n)
                            
                            for nn in neigh_neighbours:
                                if len(neigh_neighbours[nn]) < 2: continue
                                genomes = []
                                for tn in neigh_neighbours[nn]:
                                    genomes += G.node[tn]['members']
                                if len(genomes) > len(set(genomes)): continue
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
                        node_count)
        while (len(temp_c) > 0):
            G = merge_nodes(G,
                            node_count,
                            temp_c.pop(),
                            node_count + 1)
            node_count += 1

    return (G)


def clean_misassembly_edges(G, threshold):

    bad_edges = []
    for edge in G.edges():
        if float(G.edges[edge]['weight']) < (threshold * min(
                int(G.node[edge[0]]['size']), int(G.node[edge[1]]['size']))):
            bad_edges.append(edge)

    for edge in bad_edges:
        G.remove_edge(edge[0], edge[1])

    return (G)


def identify_possible_highly_variable(G,
                                   cycle_threshold_max=20,
                                   cycle_threshold_min=5,
                                   size_diff_threshold=0.5):

    # add family paralog attribute to nodes
    for node in G.nodes():
        G.node[node]['highVar'] = 0

    # find all the cycles shorter than cycle_threshold
    complete_basis = []
    for c in nx.connected_components(G):
        sub_G = G.subgraph(c)
        basis = nx.cycle_basis(sub_G, list(sub_G.nodes())[0])
        complete_basis += [
            set(b) for b in basis if len(b) <= cycle_threshold_max
        ]

    # remove cycles that are too short
    complete_basis = [
        b for b in complete_basis if len(b) >= 3
    ]

    # merge cycles with more than one node in common (nested)
    if len(complete_basis) < 1:
        return G
        
    merged_basis = [[1, set(complete_basis[0])]]
    for b in complete_basis[1:]:
        b = set(b)
        merged=False
        for i, mb in enumerate(merged_basis):
            if len(mb[1].intersection(b)) > 1:
                merged=True
                merged_basis[i][0] += 1
                merged_basis[i][1] |= b
        if not merged:
            merged_basis.append([1, b])

    for b in merged_basis:
        if b[0] < cycle_threshold_min: continue
        max_size = max([G.node[node]['size'] for node in b[1]])
        for node in b[1]:
            if G.node[node]['size'] < (size_diff_threshold*max_size):
                G.node[node]['highVar'] = 1

    return G
