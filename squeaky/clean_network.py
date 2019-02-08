import networkx as nx

# Genes at the end of contigs are more likely to be false positives thus
# we can remove those with low support
def trim_low_support_trailing_ends(G, min_support=2, max_recursive=2):

    for iter in range(max_recursive):
        bad_nodes = []
        for (node, val) in G.degree():
            if val<=1: # trailing node
                if G.node[node]['size'] < min_support:
                    bad_nodes.append(node)
        for node in bad_nodes:
            G.remove_node(node)

    return G


def identify_paralogs(G):



    return G
