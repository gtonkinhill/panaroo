import networkx as nx
from panaroo.cdhit import *
from panaroo.merge_nodes import *
from panaroo.isvalid import del_dups
from collections import defaultdict, deque, Counter
from panaroo.cdhit import is_valid
from itertools import chain, combinations
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from scipy.sparse.csgraph import connected_components, shortest_path
from scipy.stats import mode
from tqdm import tqdm
from intbitset import intbitset
import sys


# Genes at the end of contigs are more likely to be false positives thus
# we can remove those with low support
def trim_low_support_trailing_ends(G, min_support=3, max_recursive=2):

    # fix trailing
    for i in range(max_recursive):
        bad_nodes = []
        removed = False
        for (node, val) in G.degree():
            if val <= 1:  # trailing node
                if G.nodes[node]['size'] < min_support:
                    bad_nodes.append(node)
        for node in bad_nodes:
            G.remove_node(node)
            removed = True

        if not removed: break

    return G


def mod_bfs_edges(G, source, depth_limit=None):
    """Iterate over edges in a breadth-first search.
    Modified version of 'generic_bfs_edges' from networkx
    """
    neighbors = G.neighbors

    visited = {source}
    if depth_limit is None:
        depth_limit = len(G)
    queue = deque([(source, depth_limit, neighbors(source))])
    while queue:
        parent, depth_now, children = queue[0]
        try:
            child = next(children)
            if child not in visited:
                yield parent, child, depth_now
                visited.add(child)
                if depth_now > 1:
                    queue.append((child, depth_now - 1, neighbors(child)))
        except StopIteration:
            queue.popleft()


def single_linkage(G, distances_bwtn_centroids, centroid_to_index, neighbours):
    index = []
    neigh_array = []
    for neigh in neighbours:
        for sid in G.nodes[neigh]['centroid']:
            index.append(centroid_to_index[sid])
            neigh_array.append(neigh)
    index = np.array(index, dtype=int)
    neigh_array = np.array(neigh_array)

    n_components, labels = connected_components(
        csgraph=distances_bwtn_centroids[index][:, index],
        directed=False,
        return_labels=True)
    # labels = labels[index]
    for neigh in neighbours:
        l = list(set(labels[neigh_array == neigh]))
        if len(l) > 1:
            for i in l[1:]:
                labels[labels == i] = l[0]

    clusters = [
        del_dups(list(neigh_array[labels == i])) for i in np.unique(labels)
    ]

    return (clusters)


# @profile
def collapse_families(G,
                      seqid_to_centroid,
                      outdir,
                      family_threshold=0.7,
                      dna_error_threshold=0.99,
                      family_len_dif_percent=0,
                      correct_mistranslations=False,
                      length_outlier_support_proportion=0.01,
                      n_cpu=1,
                      quiet=False,
                      distances_bwtn_centroids=None,
                      centroid_to_index=None,
                      depths = [1, 2, 3],
                      search_genome_ids = None):

    node_count = max(list(G.nodes())) + 10

    if correct_mistranslations:
        threshold = [0.99, 0.98, 0.95, 0.9]
    else:
        threshold = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5]

    # precluster for speed
    if correct_mistranslations:
        cdhit_clusters = iterative_cdhit(G,
                                         outdir,
                                         thresholds=threshold,
                                         s=family_len_dif_percent,
                                         n_cpu=n_cpu,
                                         quiet=True,
                                         dna=True,
                                         word_length=7,
                                         accurate=False)
        distances_bwtn_centroids, centroid_to_index = pwdist_edlib(
            G, cdhit_clusters, dna_error_threshold, dna=True, n_cpu=n_cpu)
    elif distances_bwtn_centroids is None:
        cdhit_clusters = iterative_cdhit(G,
                                         outdir,
                                         thresholds=threshold,
                                         s=family_len_dif_percent,
                                         n_cpu=n_cpu,
                                         quiet=True,
                                         dna=False)
        distances_bwtn_centroids, centroid_to_index = pwdist_edlib(
            G, cdhit_clusters, family_threshold, dna=False, n_cpu=n_cpu)

    # keep track of centroids for each sequence. Need this to resolve clashes
    seqid_to_index = {}
    for node in G.nodes():
        for sid in G.nodes[node]['seqIDs']:
            if "refound" in sid:
                seqid_to_index[sid] = centroid_to_index[G.nodes[node]
                                                        ["longCentroidID"][1]]
            else:
                seqid_to_index[sid] = centroid_to_index[seqid_to_centroid[sid]]

    nonzero_dist = distances_bwtn_centroids.nonzero()
    nonzero_dist = set([(i, j)
                        for i, j in zip(nonzero_dist[0], nonzero_dist[1])])

    node_mem_index = {}
    for n in G.nodes():
        node_mem_index[n] = defaultdict(set)
        for sid in G.nodes[n]['seqIDs']:
            node_mem_index[n][int(sid.split("_")[0])].add(seqid_to_index[sid])

    for depth in depths:
        if not quiet: print("Processing depth: ", depth)
        if search_genome_ids is None:
            search_space = set(G.nodes())
        else:
            search_space = set()
            search_genome_ids = intbitset(search_genome_ids)
            for n in G.nodes():
                if len(G.nodes[n]['members'].intersection(search_genome_ids))>0:
                    search_space.add(n)
            
        iteration_num = 1
        while len(search_space) > 0:
            # look for nodes to merge
            temp_node_list = list(search_space)
            removed_nodes = set()
            if not quiet: print("Iteration: ", iteration_num)
            iteration_num += 1
            for node in tqdm(temp_node_list, disable=quiet):
                if node in removed_nodes: continue

                if G.degree[node] <= 2:
                    search_space.remove(node)
                    removed_nodes.add(node)
                    continue

                # find neighbouring nodes and cluster their centroid with cdhit
                neighbours = [
                    v
                    for u, v in nx.bfs_edges(G, source=node, depth_limit=depth)
                ] + [node]

                # find clusters
                clusters = single_linkage(G, distances_bwtn_centroids,
                                          centroid_to_index, neighbours)

                for cluster in clusters:

                    # check if there are any to collapse
                    if len(cluster) <= 1: continue

                    # check for conflicts
                    seen = G.nodes[cluster[0]]['members'].copy()
                    noconflict = True
                    for n in cluster[1:]:
                        if not seen.isdisjoint(G.nodes[n]['members']):
                            noconflict = False
                            break
                        seen |= G.nodes[n]['members']

                    if noconflict:
                        # no conflicts so merge
                        node_count += 1
                        for neig in cluster:
                            removed_nodes.add(neig)
                            if neig in search_space: search_space.remove(neig)

                        G = merge_node_cluster(
                            G,
                            cluster,
                            node_count,
                            multi_centroid=(not correct_mistranslations))

                        node_mem_index[node_count] = node_mem_index[cluster[0]]
                        for n in cluster[1:]:
                            for m in node_mem_index[n]:
                                node_mem_index[node_count][
                                    m] |= node_mem_index[n][m]
                            node_mem_index[n].clear()
                            node_mem_index[n] = None

                        search_space.add(node_count)
                    else:
                        # merge if the centroids don't conflict and the nodes are adjacent in the conflicting genome
                        # this corresponds to a mistranslation/frame shift/premature stop where one gene has been split
                        # into two in a subset of genomes

                        # sort by size
                        cluster = sorted(cluster,
                                         key=lambda x: G.nodes[x]['size'],
                                         reverse=True)

                        node_mem_count = Counter(
                            itertools.chain.from_iterable(
                                gen_node_iterables(G, cluster, 'members')))
                        mem_count = np.array(list(node_mem_count.values()))
                        merge_same_members = True
                        if np.sum(mem_count == 1) / float(
                                len(mem_count
                                    )) < length_outlier_support_proportion:
                            # do not merge nodes that have the same members as this is likely to be a spurious long gene
                            merge_same_members = False

                        while len(cluster) > 0:
                            sub_clust = [cluster[0]]
                            nA = cluster[0]
                            for nB in cluster[1:]:
                                mem_inter = list(
                                    G.nodes[nA]['members'].intersection(
                                        G.nodes[nB]['members']))
                                if len(mem_inter) > 0:
                                    if merge_same_members:
                                        shouldmerge = True
                                        if len(
                                                set(G.nodes[nA]['centroid']).
                                                intersection(
                                                    set(G.nodes[nB]
                                                        ['centroid']))) > 0:
                                            shouldmerge = False

                                        if shouldmerge:
                                            edge_mem_count = Counter()
                                            for e in itertools.chain.from_iterable(
                                                    gen_edge_iterables(
                                                        G, G.edges([nA, nB]),
                                                        'members')):
                                                edge_mem_count[e] += 1
                                                if edge_mem_count[e] > 3:
                                                    shouldmerge = False
                                                    break

                                        if shouldmerge:
                                            for imem in mem_inter:
                                                for sidA in node_mem_index[nA][
                                                        imem]:
                                                    for sidB in node_mem_index[
                                                            nB][imem]:
                                                        if ((
                                                                sidA, sidB
                                                        ) in nonzero_dist) or (
                                                            (sidB, sidA) in
                                                                nonzero_dist):
                                                            shouldmerge = False
                                                            break
                                                    if not shouldmerge: break
                                                if not shouldmerge: break

                                        if shouldmerge:
                                            sub_clust.append(nB)
                                else:
                                    sub_clust.append(nB)

                            if len(sub_clust) > 1:

                                clique_clusters = single_linkage(
                                    G, distances_bwtn_centroids,
                                    centroid_to_index, sub_clust)
                                for clust in clique_clusters:
                                    if len(clust) <= 1: continue
                                    node_count += 1
                                    for neig in clust:
                                        removed_nodes.add(neig)
                                        if neig in search_space:
                                            search_space.remove(neig)
                                    G = merge_node_cluster(
                                        G,
                                        clust,
                                        node_count,
                                        multi_centroid=(
                                            not correct_mistranslations),
                                        check_merge_mems=False)

                                    node_mem_index[
                                        node_count] = node_mem_index[clust[0]]
                                    for n in clust[1:]:
                                        for m in node_mem_index[n]:
                                            node_mem_index[node_count][
                                                m] |= node_mem_index[n][m]
                                        node_mem_index[n].clear()
                                        node_mem_index[n] = None

                                    search_space.add(node_count)

                            cluster = [
                                n for n in cluster if n not in sub_clust
                            ]

                if node in search_space:
                    search_space.remove(node)

    return G, distances_bwtn_centroids, centroid_to_index


def collapse_paralogs(G, centroid_contexts, max_context=5, quiet=False):

    node_count = max(list(G.nodes())) + 10

    # first sort by context length, context dist to ensure ties
    #  are broken the same way
    for centroid in centroid_contexts:
        centroid_contexts[centroid] = sorted(centroid_contexts[centroid])

    # set up for context search
    centroid_to_index = {}
    ncentroids = -1
    for node in G.nodes():
        centroid = G.nodes[node]['centroid'][0]
        if centroid not in centroid_to_index:
            ncentroids += 1
            centroid_to_index[centroid] = ncentroids
            centroid_to_index[G.nodes[node]['centroid'][0]] = ncentroids
        else:
            centroid_to_index[G.nodes[node]['centroid']
                              [0]] = centroid_to_index[centroid]
    ncentroids += 1

    for centroid in tqdm(centroid_contexts, disable=quiet):
        # calculate distance
        member_paralogs = defaultdict(list)
        for para in centroid_contexts[centroid]:
            member_paralogs[para[1]].append(para)

        ref_paralogs = max(member_paralogs.items(), key=lambda x: len(x[1]))[1]
        # for each paralog find its closest reference paralog
        cluster_dict = defaultdict(set)
        cluster_mems = defaultdict(set)
        for c, ref in enumerate(ref_paralogs):
            cluster_dict[c].add(ref[0])
            cluster_mems[c].add(ref[1])

        for para in centroid_contexts[centroid]:
            d_max = np.inf
            s_max = -np.inf
            best_cluster = None

            if para[1] == ref_paralogs[0][1]:
                # this is the reference so skip
                continue

            # first attempt by shortest path
            for c, ref in enumerate(ref_paralogs):
                if para[1] in cluster_mems[c]:
                    #dont match paralogs of the same isolate
                    continue
                # d = spath[para[0], ref[0]]
                # d = gt.shortest_distance(Gt, para[0], ref[0])
                try:
                    d = nx.shortest_path_length(G, ref[0], para[0])
                except nx.NetworkXNoPath:
                    continue
                if d < d_max:
                    d_max = d
                    best_cluster = c

            # if this fails use context
            if d_max == np.inf:
                best_cluster = 0
                s_max = -np.inf
                para_context = np.zeros(ncentroids)
                for u, node, depth in mod_bfs_edges(G, para[0], max_context):
                    para_context[centroid_to_index[G.nodes[node]['centroid']
                                                   [0]]] = depth
                for c, ref in enumerate(ref_paralogs):
                    if para[1] in cluster_mems[c]:
                        #dont match paralogs of the same isolate
                        continue
                    ref_context = np.zeros(ncentroids)
                    for u, node, depth in mod_bfs_edges(
                            G, ref[0], max_context):
                        ref_context[centroid_to_index[G.nodes[node]['centroid']
                                                      [0]]] = depth
                    s = np.sum(1 / (1 + np.abs((para_context - ref_context)[
                        (para_context * ref_context) != 0])))
                    if s > s_max:
                        s_max = s
                        best_cluster = c

            cluster_dict[best_cluster].add(para[0])
            cluster_mems[best_cluster].add(para[1])

        # merge
        for cluster in cluster_dict:
            if len(cluster_dict[cluster]) < 2: continue
            node_count += 1

            G = merge_node_cluster(G, list(cluster_dict[cluster]), node_count)

    return (G)


def merge_paralogs(G):

    node_count = max(list(G.nodes())) + 10

    # group paralog nodes by centroid
    paralog_centroids = defaultdict(list)
    for node in G.nodes():
        if G.nodes[node]['paralog']:
            for centroid in G.nodes[node]['centroid']:
                paralog_centroids[centroid].append(node)

    # find nodes that share common centroids
    paralog_centroids = paralog_centroids.values()
    merge_clusters = []
    while len(paralog_centroids) > 0:
        first, *rest = paralog_centroids
        first = set(first)
        lf = -1
        while len(first) > lf:
            lf = len(first)
            rest2 = []
            for r in rest:
                if len(first.intersection(set(r))) > 0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2
        merge_clusters.append(first)
        paralog_centroids = rest

    # merge paralog nodes that share the same centroid
    for temp_c in merge_clusters:
        if len(temp_c) > 1:
            node_count += 1
            G = merge_node_cluster(G,
                                   temp_c,
                                   node_count,
                                   check_merge_mems=False)

    return (G)


def clean_misassembly_edges(G, edge_support_threshold):

    bad_edges = set()
    max_weight = 0

    # remove edges with low support near contig ends
    for node in G.nodes():
        max_weight = max(max_weight, G.nodes[node]['size'])
        for neigh in G.neighbors(node):
            if G.nodes[neigh]['hasEnd']:
                if G[node][neigh]['size'] < edge_support_threshold:
                    bad_edges.add((node, neigh))

    # remove edges that have much lower support than the nodes they connect
    for edge in G.edges():
        if float(G.edges[edge]['size']) < (0.05 * min(
                int(G.nodes[edge[0]]['size']), int(G.nodes[edge[1]]['size']))):
            if float(G.edges[edge]['size']) < edge_support_threshold:
                bad_edges.add(edge)

    for edge in bad_edges:
        if G.has_edge(edge[0], edge[1]):
            G.remove_edge(edge[0], edge[1])

    return (G)


def identify_possible_highly_variable(G,
                                      cycle_threshold_max=20,
                                      cycle_threshold_min=5,
                                      size_diff_threshold=0.5):

    # add family paralog attribute to nodes
    for node in G.nodes():
        G.nodes[node]['highVar'] = 0

    # find all the cycles shorter than cycle_threshold
    complete_basis = []
    for c in nx.connected_components(G):
        sub_G = G.subgraph(c)
        basis = nx.cycle_basis(sub_G, list(sub_G.nodes())[0])
        complete_basis += [
            set(b) for b in basis if len(b) <= cycle_threshold_max
        ]

    # remove cycles that are too short
    complete_basis = [b for b in complete_basis if len(b) >= 3]

    # merge cycles with more than one node in common (nested)
    if len(complete_basis) < 1:
        return G

    merged_basis = [[1, set(complete_basis[0])]]
    for b in complete_basis[1:]:
        b = set(b)
        merged = False
        for i, mb in enumerate(merged_basis):
            if len(mb[1].intersection(b)) > 1:
                merged = True
                merged_basis[i][0] += 1
                merged_basis[i][1] |= b
        if not merged:
            merged_basis.append([1, b])

    for b in merged_basis:
        if b[0] < cycle_threshold_min: continue
        max_size = max([G.nodes[node]['size'] for node in b[1]])
        for node in b[1]:
            if G.nodes[node]['size'] < (size_diff_threshold * max_size):
                G.nodes[node]['highVar'] = 1

    return G
