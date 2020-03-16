import itertools
from collections import Counter
from .isvalid import del_dups
import numpy as np

def gen_node_iterables(G, nodes, feature, split=None):
    for n in nodes:
        if split is None:
            yield G.nodes[n][feature]
        else:
            yield G.nodes[n][feature].split(split)

def temp_iter(list_list):
    for n in list_list:
        yield n

def iter_del_dups(iterable):
    seen = {}
    for f in itertools.chain.from_iterable(iterable):
        seen[f] = None
    return (list(seen.keys()))

def merge_node_cluster(G,
                nodes,
                newNode,
                multi_centroid=True,
                check_merge_mems=True,
                filter_outliers=False,
                len_support_prop=0.01):

    if check_merge_mems:
        mem_count = Counter(itertools.chain.from_iterable(gen_node_iterables(G, nodes, 'members')))
        if max(mem_count.values()) > 1:
            raise ValueError("merging nodes with the same genome IDs!")


    # filter length outliers if requested
    if filter_outliers:
        lengths = np.array(list(itertools.chain.from_iterable(gen_node_iterables(G,nodes,'lengths'))))
        total_support = np.sum(list(gen_node_iterables(G,nodes,'size')))
        quantiles = np.quantile(lengths, [0.25,0.75])
        iqr15_threshold = quantiles[1] + 1.5*(quantiles[1]-quantiles[0])
        filter_nodes = []
        for node in nodes:
            if max(G.nodes[node]['lengths']) > iqr15_threshold:
                if G.nodes[node]['size'] <= len_support_prop*total_support:
                    filter_nodes.append(node)
        G.remove_nodes_from(filter_nodes)
        filter_nodes = set(filter_nodes)
        nodes = [n for n in nodes if n not in filter_nodes]

    # take node with most support as the 'consensus'
    nodes = sorted(nodes, key=lambda x: G.nodes[x]['size'])

    # First create a new node and combine the attributes
    dna = iter_del_dups(gen_node_iterables(G,nodes,'dna'))
    maxLenId = 0
    max_l = 0
    for i, s in enumerate(dna):
        if len(s) >= max_l:
            max_l = len(s)
            maxLenId = i
    members = set(iter_del_dups(gen_node_iterables(G,nodes,'members')))

    if multi_centroid:
        mergedDNA = any(gen_node_iterables(G,nodes,'mergedDNA'))
    else:
        mergedDNA = True

    G.add_node(
        newNode,
        size=len(members),
        centroid=iter_del_dups(gen_node_iterables(G,nodes,'centroid')),
        maxLenId=maxLenId,
        members=members,
        seqIDs=set(iter_del_dups(gen_node_iterables(G,nodes,'seqIDs'))),
        hasEnd=any(gen_node_iterables(G,nodes,'hasEnd')),
        protein=iter_del_dups(gen_node_iterables(G,nodes,'protein')),
        dna=dna,
        annotation=";".join(iter_del_dups(gen_node_iterables(G,nodes,'annotation',split=";"))),
        description=";".join(iter_del_dups(gen_node_iterables(G,nodes,'description',split=";"))),
        lengths=list(itertools.chain.from_iterable(gen_node_iterables(G,nodes,'lengths'))),
        longCentroidID=max(gen_node_iterables(G,nodes,'longCentroidID')),
        paralog=any(gen_node_iterables(G,nodes,'paralog')),
        mergedDNA=mergedDNA
    )
    if "prevCentroids" in G.nodes[nodes[0]]:
        G.nodes[newNode]['prevCentroids'] = ";".join(
            set(iter_del_dups(gen_node_iterables(G,nodes,'prevCentroids',split=";"))))

    # Now iterate through neighbours of each node and add them to the new node
    merge_nodes = set(nodes)
    for node in nodes:
        for neighbour in G.neighbors(node):
            if neighbour in merge_nodes: continue
            if G.has_edge(newNode, neighbour):
                G[newNode][neighbour]['weight'] += G[node][neighbour]['weight']
                G[newNode][neighbour]['members'] |= G[node][neighbour]['members']
            else:
                G.add_edge(newNode, neighbour,
                       weight=G[node][neighbour]['weight'],
                       members=G[node][neighbour]['members'])

    # remove old nodes from Graph
    G.remove_nodes_from(nodes)

    return G



def delete_node(G, node):
    # add in new edges
    for mem in G.nodes[node]['members']:
        mem_edges = list(
            set([e[1] for e in G.edges(node) if mem in G.edges[e]['members']]))
        if len(mem_edges) < 2: continue
        for n1, n2 in itertools.combinations(mem_edges, 2):
            if G.has_edge(n1, n2):
                G[n1][n2]['members'] |= set([mem])
                G[n1][n2]['weight'] += 1
            else:
                G.add_edge(n1, n2, weight=1, members=set([mem]))

    # now remove node
    G.remove_node(node)

    return G


def remove_member_from_node(G, node, member):

    # add in replacement edges if required
    mem_edges = list(
        set([e[1] for e in G.edges(node) if member in G.edges[e]['members']]))
    if len(mem_edges) > 1:
        for n1, n2 in itertools.combinations(mem_edges, 2):
            if G.has_edge(n1, n2):
                G[n1][n2]['members'] |= [member]
                G[n1][n2]['weight'] += 1
            else:
                G.add_edge(n1, n2, weight=1, members=set([member]))

    # remove member from node
    while str(member) in G.nodes[node]['members']:
        G.nodes[node]['members'].remove(str(member))
    G.nodes[node]['seqIDs'] = set([
        sid for sid in G.nodes[node]['seqIDs']
        if sid.split("_")[0] != str(member)
    ])
    G.nodes[node]['size'] -= 1

    return G