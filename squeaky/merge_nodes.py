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
