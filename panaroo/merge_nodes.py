import itertools
from panaroo.isvalid import del_dups

def merge_nodes(G,
                nodeA,
                nodeB,
                newNode,
                multi_centroid=True,
                check_merge_mems=True):

    if check_merge_mems:
        if len(set(G.nodes[nodeA]['members'])
               & set(G.nodes[nodeB]['members'])) > 0:
            raise ValueError("merging nodes with the same genome IDs!")

    # take node with most support as the 'consensus'
    if G.nodes[nodeA]['size'] < G.nodes[nodeB]['size']:
        nodeB, nodeA = nodeA, nodeB

    # First create a new node and combine the attributes
    if multi_centroid:
        G.add_node(newNode,
                   size=len(set(G.nodes[nodeA]['members'] + G.nodes[nodeB]['members'])),
                   centroid=";".join(
                       del_dups(G.nodes[nodeA]['centroid'].split(";") +
                           G.nodes[nodeB]['centroid'].split(";"))),
                   members=G.nodes[nodeA]['members'] + G.nodes[nodeB]['members'],
                   seqIDs=G.nodes[nodeA]['seqIDs'] + G.nodes[nodeB]['seqIDs'],
                   hasEnd=(G.nodes[nodeA]['hasEnd'] or G.nodes[nodeB]['hasEnd']),
                   protein=";".join(
                       del_dups(G.nodes[nodeA]['protein'].split(";") +
                           G.nodes[nodeB]['protein'].split(";"))),
                   dna=";".join(G.nodes[nodeA]['dna'].split(";") +
                                G.nodes[nodeB]['dna'].split(";")),
                   annotation=";".join(
                       del_dups(G.nodes[nodeA]['annotation'].split(";") +
                           G.nodes[nodeB]['annotation'].split(";"))),
                   description=";".join(
                       del_dups(G.nodes[nodeA]['description'].split(";") +
                           G.nodes[nodeB]['description'].split(";"))),
                   lengths=G.nodes[nodeA]['lengths'] + G.nodes[nodeB]['lengths'],
                   longCentroidID=max(G.nodes[nodeA]['longCentroidID'], G.nodes[nodeB]['longCentroidID']),
                   paralog=(G.nodes[nodeA]['paralog']
                            or G.nodes[nodeB]['paralog']),
                   mergedDNA=(G.nodes[nodeA]['mergedDNA']
                              or G.nodes[nodeB]['mergedDNA']))
        if "prevCentroids" in G.nodes[nodeA]:
            G.nodes[newNode]['prevCentroids'] = ";".join(
                       set(G.nodes[nodeA]['prevCentroids'].split(";") +
                           G.nodes[nodeB]['prevCentroids'].split(";")))
    else:
        G.add_node(newNode,
                   size=len(set(G.nodes[nodeA]['members'] + G.nodes[nodeB]['members'])),
                   centroid=";".join(
                       del_dups(G.nodes[nodeA]['centroid'].split(";") +
                           G.nodes[nodeB]['centroid'].split(";"))),
                   members=G.nodes[nodeA]['members'] + G.nodes[nodeB]['members'],
                   seqIDs=G.nodes[nodeA]['seqIDs'] + G.nodes[nodeB]['seqIDs'],
                   hasEnd=(G.nodes[nodeA]['hasEnd'] or G.nodes[nodeB]['hasEnd']),
                   protein=";".join(
                       del_dups(G.nodes[nodeA]['protein'].split(";"))+
                                G.nodes[nodeB]['protein'].split(";")), 
                   dna=";".join(G.nodes[nodeA]['dna'].split(";") +
                                G.nodes[nodeB]['dna'].split(";")),
                   annotation=G.nodes[nodeA]['annotation'],
                   description=G.nodes[nodeA]['description'],
                   paralog=(G.nodes[nodeA]['paralog']
                            or G.nodes[nodeB]['paralog']),
                   lengths=G.nodes[nodeA]['lengths'] + G.nodes[nodeB]['lengths'],
                   longCentroidID=max(G.nodes[nodeA]['longCentroidID'], G.nodes[nodeB]['longCentroidID']),
                   mergedDNA=True)
        if "prevCentroids" in G.nodes[nodeA]:
            G.nodes[newNode]['prevCentroids'] = ";".join(
                       set(G.nodes[nodeA]['prevCentroids'].split(";") +
                           G.nodes[nodeB]['prevCentroids'].split(";")))

    # Now iterate through neighbours of each node and add them to the new node
    neigboursB = list(G.neighbors(nodeB))
    neigboursA = list(G.neighbors(nodeA))
    for neighbor in neigboursA:
        if neighbor in neigboursB:
            G.add_edge(newNode,
                       neighbor,
                       weight=G[nodeA][neighbor]['weight'] +
                       G[nodeB][neighbor]['weight'],
                       members=G[nodeA][neighbor]['members'] +
                       G[nodeB][neighbor]['members'])
            neigboursB.remove(neighbor)
        else:
            G.add_edge(newNode,
                       neighbor,
                       weight=G[nodeA][neighbor]['weight'],
                       members=G[nodeA][neighbor]['members'])

    for neighbor in neigboursB:
        G.add_edge(newNode,
                   neighbor,
                   weight=G[nodeB][neighbor]['weight'],
                   members=G[nodeB][neighbor]['members'])

    # remove old nodes from Graph
    G.remove_nodes_from([nodeA, nodeB])

    if len(max(G.nodess[newNode]["dna"].split(";"), key=len)) <= 0:
        print(G.nodes[newNode]["dna"])
        raise NameError("Problem!")

    return G


def delete_node(G, node):
    # add in new edges
    for mem in G.nodes[node]['members']:
        mem_edges = list(
            set([e[1] for e in G.edges(node) if mem in G.edges[e]['members']]))
        if len(mem_edges)<2: continue
        for n1, n2 in itertools.combinations(mem_edges, 2):
            if G.has_edge(n1, n2):
                G[n1][n2]['members'] += [mem]
                G[n1][n2]['weight'] += 1
            else:
                G.add_edge(n1, n2, weight=1, members=[mem])

    # now remove node
    G.remove_node(node)

    return G


def remove_member_from_node(G, node, member):

    # add in replacement edges if required
    mem_edges = list(
            set([e[1] for e in G.edges(node) if member in G.edges[e]['members']]))
    if len(mem_edges)>1:
        for n1, n2 in itertools.combinations(mem_edges, 2):
            if G.has_edge(n1, n2):
                G[n1][n2]['members'] += [member]
                G[n1][n2]['weight'] += 1
            else:
                G.add_edge(n1, n2, weight=1, members=[member])

    # remove member from node
    while str(member) in G.nodes[node]['members']:
        G.nodes[node]['members'].remove(str(member))
    G.nodes[node]['seqIDs'] = [sid for sid in G.nodes[node]['seqIDs'] if sid.split("_")[0]!=str(member)]
    G.nodes[node]['size'] -= 1

    return G