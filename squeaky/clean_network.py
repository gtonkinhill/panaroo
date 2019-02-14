import networkx as nx

# Genes at the end of contigs are more likely to be false positives thus
# we can remove those with low support
def trim_low_support_trailing_ends(G, min_support=3, max_recursive=2):

    # fix trailing
    for iter in range(max_recursive):
        bad_nodes = []
        for (node, val) in G.degree():
            if val<=1: # trailing node
                if G.node[node]['size'] < min_support:
                    bad_nodes.append(node)
        for node in bad_nodes:
            G.remove_node(node)

        # fix trailing due to paralog
        bad_nodes = []
        for node in G.nodes():
            centroids = []
            for neighbor in G.neighbors(node):
                centroids.append(G.node[neighbor]['centroid'])
            if len(set(centroids))<=1:
                if len(centroids) < min_support:
                    bad_nodes.append(node)
        for node in bad_nodes:
            G.remove_node(node)

    return G


# Look for cycles and collapse nodes within those that are within
# family_threshold pairwise sequence identity
def collapse_families(G, cycle_threshold, family_threshold, outdir,
    dna_error_threshold=0.95,
    correct_mistranslations=True):

    node_count = max(G.nodes())+1

    # find all the cycles shorter than cycle_threshold
    complete_basis = []
    for c in nx.connected_components(G):
        sub_G = G.subgraph(c)
        basis = nx.cycle_basis(sub_G, sub_G.nodes()[0])
        complete_basis += [set(b) for b in basis if len(b) < cycle_threshold]

    # remove paralogs and resulting cycles that are too short
    temp_basis = []
    for i, b in enumerate(complete_basis):
        new_b = set()
        for node in b:
            if not G.node[node]["paralog"]:
                new_b.add(node)
        if len(new_b)>1:
            temp_basis.append(new_b)
    complete_basis = temp_basis

    # merge cycles with more than one node in common (nested)
    merged_basis = []
    while len(complete_basis)>0:
        first, *rest = complete_basis
        is_merged = False
        while not is_merged:
            is_merged = True
            rest2 = []
            for r in rest:
                if len(first.intersection(r)>1):
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
        for b in merged_basis:
            clusters = cluster_nodes_cdhit(G, b, outdir, id=dna_error_threshold,
                quiet=True, dna=True)
            # now merge nodes that clustered
            temp_b = set()
            for c in clusters:
                if len(c)>1:
                    # keep the centroid with the highest support
                    top_size = -1
                    for node in c:
                        if G.node[node]["size"]>top_size:
                            top=node
                            top_size=G.node[node]["size"]
                    temp_c = c
                    if top==c[0]:
                        G = merge_nodes(G, top, c[1], node_count)
                        temp_c.remove(top)
                        temp_c.remove(c[1])
                    else:
                        G = merge_nodes(G, top, c[0], node_count)
                        temp_c.remove(top)
                        temp_c.remove(c[0])
                    while (len(temp_c)>0):
                        G = merge_nodes(G, node_count, temp_c[0], node_count+1)
                        node_count+=1
                        temp_c.remove(temp_c[0])
                    temp_b.add(node_count-1)
                else:
                    temp_b.add(c)
            cleaned_merged_basis.append(temp_b)
    else:
        cleaned_merged_basis = merged_basis

    # merge nodes based on the family_threshold by clustering at the protein
    # level
    for b in merged_basis:
        clusters = cluster_nodes_cdhit(G, b, outdir, id=family_threshold,
            quiet=True, dna=False)
        # now merge nodes that clustered
        for c in clusters:
            if len(c)>1:
                # keep the centroid with the highest support
                temp_c = c
                G = merge_nodes(G, c[0], c[1], node_count,
                        multi_centroid=True)
                temp_c.remove(c[0])
                temp_c.remove(c[1])
                while (len(temp_c)>0):
                    G = merge_nodes(G, node_count, temp_c[0], node_count+1)
                    node_count+=1


    return G
