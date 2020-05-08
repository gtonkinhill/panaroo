import networkx as nx
from collections import deque, defaultdict

def conv_list(maybe_list):
    if not isinstance(maybe_list, list):
        maybe_list = [maybe_list]
    return (maybe_list)

def get_target(g, gene):
    for n, attr in g.nodes(data=True):
        if attr["name"] == gene:
            return n
    raise NameError("Gene ID does not match any in the graph!")

def bfs_with_dist(G, source, depth_limit=None, prev_depth=0, genome=None):
    successors = G.neighbors
    for e in generic_bfs_edges_with_dist(G, source, successors, depth_limit, genome=genome):
        yield e

def generic_bfs_edges_with_dist(G, source, neighbors=None, depth_limit=None, genome=None):
    visited = {source}
    if depth_limit is None:
        depth_limit = len(G)
    queue = deque([(source, depth_limit, get_neighbours_with_genome(G, source, genome))])
    while queue:
        parent, depth_now, children = queue[0]
        try:
            child = next(children)
            if child not in visited:
                yield parent, child, depth_limit-depth_now+1
                visited.add(child)
                if depth_now > 1:
                    queue.append((child, depth_now - 1, get_neighbours_with_genome(G, child, genome)))
        except StopIteration:
            queue.popleft()

def get_neighbours_with_genome(G, node, genome):
    neighbours = G.neighbors(node)
    if genome is None:
        return (neighbours)
    else:
        return (n for n in neighbours if genome in conv_list(G[node][n]['members']))

def get_options():
    import argparse

    description = 'Explore gene neighbourhood'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo-gene-neighbourhood')

    parser.add_argument("--gene",
                        type=str,
                        required=True,
                        help="gene of interest")
    parser.add_argument("--genome_id",
                        type=str,
                        default=None,
                        help="genome ID of interest (default=ALL)")
    parser.add_argument("--graph",
                        type=str,
                        required=True,
                        help="genome graph gml ('final_graph.gml')")
    parser.add_argument("--expand_no",
                        default=5,
                        help=("lengths of the path that will be expanded on" +
                              " in a radius the target gene (default=5)"),
                        type=int)
    parser.add_argument("--out", help="output file")

    args = parser.parse_args()

    return (args)

def main():
    args = get_options()

    # load graph
    G = nx.read_gml(args.graph)
    for n in G.nodes():
        G.nodes[n]['members'] = set(conv_list(G.nodes[n]['members']))

    # find target gene
    target = get_target(G, args.gene)

    # find target genome id if requested
    gid = None
    if args.genome_id is not None:
        for i, genome in enumerate(G.graph['isolateNames']):
            if genome==args.genome_id:
                gid = str(i)
        if gid is None:
            raise NameError("Genome ID does not match any in the graph!")

    # write out neighbouring genes and distance from target
    
    # allocate edges to members
    mems_to_edges = defaultdict(list)
    if gid is None:
        msearch = G.nodes[target]['members']
    else:
        msearch = [gid]

    for mem in msearch:
        for u,v,d in bfs_with_dist(G, target, depth_limit=args.expand_no, genome=mem): 
            mems_to_edges[mem].append((u,v))

    # find path for each member
    paths_to_members = defaultdict(list)
    for mem in mems_to_edges:
        # create temporary graph
        tG = nx.Graph()
        tG.add_edges_from(mems_to_edges[mem])

        # find largest connected component that contains the target
        for c in sorted(nx.connected_components(tG), key=len, reverse=True):
            if target in c:
                path = sorted(c)
                break
        
        # reorder path
        for n in path:
            if tG.degree(n)==1:
                path = [n] + [v for u, v in nx.bfs_edges(tG, source=n)]
                break
        
        paths_to_members[tuple(path)].append(mem)

    with open(args.out, 'w') as outfile:
        outfile.write("members\tpath\n")
        for path in paths_to_members:
            outfile.write(",".join([G.graph['isolateNames'][m] for m in paths_to_members[path]]) + "\t")
            outfile.write(",".join([G.nodes[n]['name'] for n in path]) + "\n")

    return


if __name__ == '__main__':
    main()
