import itertools


def merge_nodes(G,
                nodeA,
                nodeB,
                newNode,
                multi_centroid=True,
                check_merge_mems=True):

    if check_merge_mems:
        if len(set(G.node[nodeA]['members'])
               & set(G.node[nodeB]['members'])) > 0:
            raise ValueError("merging nodes with the same genome IDs!")

    # First create a new node and combine the attributes
    if multi_centroid:
        G.add_node(newNode,
                   size=G.node[nodeA]['size'] + G.node[nodeB]['size'],
                   centroid=";".join(
                       set(G.node[nodeA]['centroid'].split(";") +
                           G.node[nodeB]['centroid'].split(";"))),
                   members=G.node[nodeA]['members'] + G.node[nodeB]['members'],
                   seqIDs=G.node[nodeA]['seqIDs'] + G.node[nodeB]['seqIDs'],
                   hasEnd=(G.node[nodeA]['hasEnd'] or G.node[nodeB]['hasEnd']),
                   protein=";".join(
                       set(G.node[nodeA]['protein'].split(";") +
                           G.node[nodeB]['protein'].split(";"))),
                   dna=";".join(G.node[nodeA]['dna'].split(";") +
                                G.node[nodeB]['dna'].split(";")),
                   annotation=";".join(
                       set(G.node[nodeA]['annotation'].split(";") +
                           G.node[nodeB]['annotation'].split(";"))),
                   description=";".join(
                       set(G.node[nodeA]['description'].split(";") +
                           G.node[nodeB]['description'].split(";"))),
                   lengths=G.node[nodeA]['lengths'] + G.node[nodeB]['lengths'],
                   longCentroidID=max(G.node[nodeA]['longCentroidID'], G.node[nodeB]['longCentroidID']),
                   paralog=(G.node[nodeA]['paralog']
                            or G.node[nodeB]['paralog']),
                   mergedDNA=(G.node[nodeA]['mergedDNA']
                              or G.node[nodeB]['mergedDNA']))
        if "prevCentroids" in G.node[nodeA]:
            G.node[newNode]['prevCentroids'] = ";".join(
                       set(G.node[nodeA]['prevCentroids'].split(";") +
                           G.node[nodeB]['prevCentroids'].split(";")))
    else:
        # take node with most support as the 'consensus'
        if G.node[nodeA]['size'] < G.node[nodeB]['size']:
            nodeB, nodeA = nodeA, nodeB

        G.add_node(newNode,
                   size=len(set(G.node[nodeA]['members'] + G.node[nodeB]['members'])),
                   centroid=G.node[nodeA]['centroid'],
                   members=G.node[nodeA]['members'] + G.node[nodeB]['members'],
                   seqIDs=G.node[nodeA]['seqIDs'] + G.node[nodeB]['seqIDs'],
                   hasEnd=(G.node[nodeA]['hasEnd'] or G.node[nodeB]['hasEnd']),
                   protein=G.node[nodeA]['protein'],
                   dna=";".join(G.node[nodeA]['dna'].split(";") +
                                G.node[nodeB]['dna'].split(";")),
                   annotation=G.node[nodeA]['annotation'],
                   description=G.node[nodeA]['description'],
                   paralog=(G.node[nodeA]['paralog']
                            or G.node[nodeB]['paralog']),
                   lengths=G.node[nodeA]['lengths'] + G.node[nodeB]['lengths'],
                   longCentroidID=max(G.node[nodeA]['longCentroidID'], G.node[nodeB]['longCentroidID']),
                   mergedDNA=True)
        if "prevCentroids" in G.node[nodeA]:
            G.node[newNode]['prevCentroids'] = ";".join(
                       set(G.node[nodeA]['prevCentroids'].split(";") +
                           G.node[nodeB]['prevCentroids'].split(";")))

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

    if len(max(G.nodes[newNode]["dna"].split(";"), key=len)) <= 0:
        print(G.node[newNode]["dna"])
        raise NameError("Problem!")

    return G


def delete_node(G, node):
    # add in new edges
    for mem in G.node[node]['members']:
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
    while str(member) in G.node[node]['members']:
        G.node[node]['members'].remove(str(member))
    G.node[node]['seqIDs'] = [sid for sid in G.node[node]['seqIDs'] if sid.split("_")[0]!=str(member)]
    G.node[node]['size'] -= 1

    return G