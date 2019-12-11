import networkx as nx


def index_neighbors(g):
    g2g2n = {}
    for n, attr in G.nodess(data=True):
        genomes = attr["genomeIDs"].split(";")
        genome_dict = dict([(genome, True) for genome in genomes])
        for neighbor in g.neighbors(n):
            edges = g.edges[n, neighbor]["members"]
            if type(edges) != list:
                edges = [edges]
            for genome in edges:
                if genome in genome_dict:
                    if n not in g2g2n:
                        g2g2n[n] = {}
                    if genome not in g2g2n[n]:
                        g2g2n[n][genome] = [neighbor]
                    else:
                        g2g2n[n][genome].append(neighbor)
    return g2g2n


def expand_path(g2g2n, paths, target):
    for genome in paths:
        path = paths[genome]
        head = path[-1]
        if genome not in g2g2n[head]:
            continue
        neighbors = g2g2n[head][genome]
        if len(neighbors) == 2:
            if len(path) == 1:
                if neighbors[0] != target:
                    path.append(neighbors[0])
                else:
                    path.append(neighbors[1])
            elif neighbors[0] != path[-2]:
                path.append(neighbors[0])
            else:
                path.append(neighbors[1])
    return paths


def join_paths(path_f, path_r, target):
    all_paths = []
    path_dict = {}
    for g in sorted(
            list(set([i for i in path_f.keys()] + [i
                                                   for i in path_r.keys()]))):
        p = []
        if g in path_f and g in path_r:
            p = path_r[g][::-1] + [target] + path_f[g]
        elif g in path_f:
            p = [target] + path_f[g]
        else:
            p = [target] + path_r[g]
        if not tuple(p) in path_dict:
            path_dict[tuple(p)] = [g]
        else:
            path_dict[tuple(p)].append(g)
    return path_dict


def get_neighbors(name, g2g2n, graph, expand_no):
    target = g2g2n[name]
    #forward paths
    paths_f = {}
    #reverse paths
    paths_r = {}
    #iterate over all genomes that are in the target node
    for genome in target:
        #forward
        if genome not in paths_f:
            paths_f[genome] = [target[genome][0]]
        else:
            paths_f[genome].append(target[genome][0])
        #reverse (always shorter or equal size)
        if len(target[genome]) > 1:
            if genome not in paths_r:
                paths_r[genome] = [target[genome][1]]
            else:
                paths_r[genome].append(target[genome][1])
    #expand paths to both sides
    for i in range(expand_no):
        expand_path(g2g2n, paths_r, name)
        expand_path(g2g2n, paths_f, name)
    #aggregate paths
    path_dict = join_paths(paths_f, paths_r, name)
    path_dict_t = {}
    #translate gene ids into gene names
    for p, genomes in zip(path_dict.keys(), path_dict.values()):
        gene_path = []
        for gene in p:
            gene_path.append(graph.nodes[gene]["name"])
        path_dict_t[tuple(genomes)] = gene_path
    return path_dict_t


def get_target(g, gene):
    for n, attr in G.nodess(data=True):
        if attr["name"] == gene:
            return n


def write_paths(paths, out):
    with open(out, 'w') as f:
        for genomes, p in zip(paths.keys(), paths.values()):
            f.write("\t".join(p))
            f.write('\n')
            f.write("\t".join(genomes))
            f.write('\n')


def run(gene, graph, expand_no, out):
    g = nx.read_gml(graph)
    target = get_target(g, gene)
    g2g2n = index_neighbors(g)
    paths = get_neighbors(target, g2g2n, g, expand_no)
    write_paths(paths, out)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser("explore gene neighbordhood")
    parser.add_argument("gene", help="gene of interested")
    parser.add_argument("graph", help="genome graph gml")
    parser.add_argument("--expand_no", default = 5, help="lengths of the path that will be expanded on either side of the target gene", type = int)
    parser.add_argument("out", help="output file")
    args = parser.parse_args()
    run(**vars(args))
