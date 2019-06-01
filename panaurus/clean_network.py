import networkx as nx
from panaurus.cdhit import cluster_nodes_cdhit
from panaurus.merge_nodes import merge_nodes
from collections import defaultdict
from panaurus.cdhit import is_valid


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


# Look for cycles and collapse nodes within those that are within
# family_threshold pairwise sequence identity
def collapse_families(G,
                      cycle_threshold,
                      family_threshold,
                      outdir,
                      dna_error_threshold=0.99,
                      correct_mistranslations=True,
                      quiet=False):

    node_count = max(list(G.nodes())) + 10

    # find all the cycles shorter than cycle_threshold
    complete_basis = []
    for c in nx.connected_components(G):
        sub_G = G.subgraph(c)
        basis = nx.cycle_basis(sub_G, list(sub_G.nodes())[0])
        complete_basis += [set(b) for b in basis if len(b) < cycle_threshold]

    # remove cycles that are too short
    complete_basis = [b for b in complete_basis if len(b) > 1]

    # merge cycles with more than one node in common (nested)
    merged_basis = []
    while len(complete_basis) > 0:
        first, *rest = complete_basis
        is_merged = False
        while not is_merged:
            is_merged = True
            rest2 = []
            for r in rest:
                if len(first.intersection(r)) > 1:
                    first |= set(r)
                    is_merged = False
                else:
                    rest2.append(r)
            rest = rest2
        merged_basis.append(first)
        complete_basis = rest

    # correct possible mistranslations using clustering at the DNA level if
    # requested.
    cleaned_merged_basis = []
    if correct_mistranslations:
        while len(merged_basis) > 0:
            b = merged_basis.pop()
            # print("basis: ", b)
            clusters = cluster_nodes_cdhit(G,
                                           b,
                                           outdir,
                                           id=dna_error_threshold,
                                           quiet=True,
                                           dna=True,
                                           use_local=True,
                                           aS=0.9)
            # if not quiet:
            #     print("cycle:", b)
            #     print("clusters:", clusters)
            # now merge nodes that clustered
            temp_b = []
            swap_dict = {}
            for c in clusters:
                if len(c) > 1:
                    # keep the centroid with the highest support
                    top_size = -1
                    for node in c:
                        if G.node[node]["size"] > top_size:
                            top = node
                            top_size = G.node[node]["size"]
                    temp_c = c.copy()
                    if top == c[0]:
                        G = merge_nodes(G, top, c[1], node_count)
                        temp_c.remove(top)
                        temp_c.remove(c[1])
                    else:
                        G = merge_nodes(G, top, c[0], node_count)
                        temp_c.remove(top)
                        temp_c.remove(c[0])
                    while (len(temp_c) > 0):
                        G = merge_nodes(G, node_count, temp_c[0],
                                        node_count + 1)
                        node_count += 1
                        temp_c.remove(temp_c[0])
                    temp_b.append(node_count)
                    node_count += 1
                else:
                    temp_b.append(c[0])
            cleaned_merged_basis.append(set(temp_b))
            # update merged_basis to use new node names
            for i, b2 in enumerate(merged_basis):
                for j, c in enumerate(clusters):
                    if len(b2.intersection(c)) > 0:
                        # intersection between basis
                        for id in c:
                            if id in merged_basis[i]:
                                merged_basis[i].remove(id)
                        merged_basis[i].add(temp_b[j])
            # update cleaned_merged_basis to use new node file_names
            for i, b2 in enumerate(cleaned_merged_basis):
                for j, c in enumerate(clusters):
                    if len(b2.intersection(c)) > 0:
                        # intersection between basis
                        for id in c:
                            if id in cleaned_merged_basis[i]:
                                cleaned_merged_basis[i].remove(id)
                        cleaned_merged_basis[i].add(temp_b[j])
    else:
        cleaned_merged_basis = merged_basis

    # merge nodes based on the family_threshold by clustering at the protein
    # level
    node_count += 1
    while len(cleaned_merged_basis) > 0:
        b = cleaned_merged_basis.pop()
        # print("cleaned basis: ", b)
        clusters = cluster_nodes_cdhit(G,
                                       b,
                                       outdir,
                                       id=family_threshold,
                                       quiet=True,
                                       dna=False)
        # if not quiet:
        #     print("cycle:", b)
        #     print("fam clusters:", clusters)
        # now merge nodes that clustered
        temp_b = []
        for c in clusters:
            if len(c) > 1:
                # keep the centroid with the highest support
                temp_c = c.copy()
                G = merge_nodes(G,
                                temp_c.pop(),
                                temp_c.pop(),
                                node_count,
                                multi_centroid=True)
                while (len(temp_c) > 0):
                    G = merge_nodes(G,
                                    node_count,
                                    temp_c.pop(),
                                    node_count + 1,
                                    multi_centroid=True)
                    node_count += 1
                temp_b.append(node_count)
            else:
                temp_b.append(c[0])
            node_count += 1
        # update cleaned_merged_basis to use new node file_names
        for i, b2 in enumerate(cleaned_merged_basis):
            for j, c in enumerate(clusters):
                if len(b2.intersection(c)) > 0:
                    # intersection between basis
                    for id in c:
                        if id in cleaned_merged_basis[i]:
                            cleaned_merged_basis[i].remove(id)
                    cleaned_merged_basis[i].add(temp_b[j])

    # merge paralogs that now have the exact same context (using node number)
    # first find those with the same context
    para_clusters = defaultdict(list)
    for node, data in G.nodes(data=True):
        if data['paralog']:
            context = tuple([data['centroid']] +
                            sorted([n for n in G.neighbors(node)]))
            para_clusters[context].append(node)
    # now merge
    for context in para_clusters:
        if len(para_clusters[context]) > 1:
            G = merge_nodes(G, para_clusters[context].pop(),
                            para_clusters[context].pop(), node_count + 1)
            while len(para_clusters[context]) > 0:
                node_count += 1
                G = merge_nodes(G, para_clusters[context].pop(), node_count,
                                node_count + 1)

    return G


def collapse_paralogs(G, quiet=False):

    node_count = max(list(G.nodes())) + 10
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
                    v for u, v in nx.bfs_edges(G, source=node, depth_limit=3)
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
