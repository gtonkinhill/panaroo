import networkx as nx


def trim_low_support_trailing_ends(G, min_support=2, max_recursive=2):

    for iter in range(max_recursive):
        for (node, val) in G.degree():
            if val<=1: # trailing node
                if G.node[node]['size'] < min_support:
                    G.remove_node(node)

    return G
