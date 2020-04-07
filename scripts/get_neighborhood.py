import networkx as nx
from collections import deque

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
                yield child, depth_limit-depth_now+1
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
        return (n for n in neighbours if genome in G[node][n]['members'])


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
        G.nodes[n]['members'] = set(G.nodes[n]['members'])

    # find target gene
    target = get_target(G, args.gene)

    # find target genome id if requested
    if args.genome_id is not None:
        gid = None
        for i, genome in enumerate(G.graph['isolateNames']):
            if genome==args.genome_id:
                gid = str(i)
        if gid is None:
            raise NameError("Genome ID does not match any in the graph!")

    # write out neighbouring genes and distance from target
    with open(args.out, 'w') as outfile:
        outfile.write("gene,dist\n")
        for node, dist in bfs_with_dist(G, target, depth_limit=args.expand_no, genome=gid):
            gname = G.nodes[node]['name']
            outfile.write(gname + "," + str(dist) +"\n")

    return


if __name__ == '__main__':
    main()
