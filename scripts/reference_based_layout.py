import networkx as nxvisited
import re
import networkx.classes.function as function
import networkx.algorithms.connectivity.cuts as cuts
import networkx as nx
import pandas as pd


def get_dist(ref_s, ref_t, max_dist):
    s = ref_s.split("_")[-1]
    t = ref_t.split("_")[-1]
    return min(abs(int(s) - int(t)), abs(abs(int(s) - int(t)) - max_dist))


def add_to_queue(g, s, nodes, visited, sink, mapping, ref_g_id, max_dist):
    add = []
    for i in nodes:
        if i in visited:
            continue
        #if "101_0_" not in g.nodes[i]['attr_dict']['geneIDs']: #DCC 1/2
        #if "1419_0_" not in g.nodes[i]['attr_dict']['geneIDs']: # K3
        #if "5_0_" not in g.nodes[i]['attr_dict']['geneIDs']: #tb relatives
        if "%s_0_" % ref_g_id not in g.nodes[i]['geneIDs'] and i not in visited:
            add.append(i)
        else:
            #print(g.nodes[s]['attr_dict']['name'],g.nodes[i]['attr_dict']['name'] )
            dist = get_dist(mapping.loc[g.nodes[s]['name'], "gene_id"],
                            mapping.loc[g.nodes[i]['name'], "gene_id"],
                            max_dist)
            #if sink["sink"] is not None:
            if dist > 100:
                #print("distance is", dist)
                sink["sink"] = i
        visited.add(i)
    return add


def create_mapping(g, ref_g_id):
    #look up table for name vs node id
    gene_dict = {}
    for n in list(g):
        gene_ids = [(i.split("_")[0], i)
                    for i in g.nodes[n]['geneIDs'].split(";")]
        gene_ids = list(filter(lambda x: ref_g_id == x[0], gene_ids))
        #g.nodes[n]['attr_dict']['name'], g.nodes[n]['attr_dict']['id']) for n in list(g)])
        #e = re.compile("1419_0_[0-9]*") #DCC1/2
        #e = re.compile("101_0_[0-9]*")
        #e = re.compile("5_0_[0-9]*")
        #e = re.compile("%s_0_[0-9]*" % ref_g_id)
        #m = e.search(gene_ids)
        #if m:
        #    gene_dict[g.nodes[n]['name']] =  m.group()
        if len(gene_ids) != 0:
            gene_dict[g.nodes[n]['name']] = gene_ids[0][1]
    mapping = pd.DataFrame.from_dict(gene_dict, orient='index')
    mapping.columns = ["gene_id"]
    #we dont want to include refound genes in this step TODO also consider in add reference edges step
    mapping = mapping.loc[~mapping.loc[:, "gene_id"].str.contains("refound"), ]
    return mapping


def add_ref_edges(g, mapping):
    name_dict = dict([(g.nodes[n]['name'], n) for n in list(g)])
    for n in mapping.index:
        mapping.loc[n, "seq"] = int(mapping.loc[n, "gene_id"].split("_")[2])
    mapping.sort_values("seq", inplace=True)
    j = 0
    for i in range(1, mapping.shape[0]):
        node1 = str(name_dict[mapping.index[i]])
        node2 = str(name_dict[mapping.index[i - 1]])
        if not g.has_edge(node1, node2):
            j += 1
            g.add_edge(node1, node2)
    return g


def remove_var_edges(g):
    for n in g:
        if g.nodes[n]["highVar"] == 1 and "%s_0_" % ref_g_id not in g.nodes[n][
                'geneIDs']:
            var_nodes.append(n)
    var_nodes = []
    g.remove_nodes_from(var_nodes)
    return g


def layout(graph, ref_g_id, cut_edges_out, ignore_high_var,
           add_reference_edges):
    g = nx.read_gml(graph)
    #look up table for name vs node id
    mapping = create_mapping(g, ref_g_id)
    gene_order = [
        int(mapping.loc[n, "gene_id"].split("_")[2]) for n in mapping.index
    ]
    max_dist = max(gene_order)
    if ignore_high_var:
        g = remove_var_edges(g)
    if add_reference_edges:
        g = add_ref_edges(g, mapping)
        #write gml with reference edges to disk to be used in cytoscape instead of the original final_graph.gml
        nx.write_gml(g, "with_ref_" + graph)
    name_dict = dict([(g.nodes[n]['name'], n) for n in list(g)])
    #set capacity for edges for the min cut algorithm as the weight of that edge
    for e in g.edges:
        try:
            g.edges[e]["capacity"] = g.edges[e]["weight"]
        except:
            g.edges[e]["capacity"] = 1
    #store edges to be taken out of the graph
    cut_edges = []
    i = 0
    cur_try = 0
    #iterate over all reference nodes in mapping table
    while i < len(mapping.index):
        n = mapping.index[i]
        print(i)
        if n not in name_dict:
            i += 1
            continue
        nid = name_dict[n]
        visited = set([nid])
        sink = {"sink": None}
        queue = add_to_queue(g, nid, g.neighbors(nid), visited, sink, mapping,
                             ref_g_id, max_dist)
        #depth first search
        last_target = None
        while len(queue) != 0:
            #print(len(queue))
            target = queue.pop(0)
            visited.add(target)
            #print(len(visited))
            neighbors = g.neighbors(target)
            #for each reference node explore all edges that lead to non-reference nodes
            queue = queue + add_to_queue(g, nid, neighbors, visited, sink,
                                         mapping, ref_g_id, max_dist)
        last_target = None
        #did we find a long-range connection?
        if sink["sink"] is not None:
            print("found path")
            visited.add(sink["sink"])
            s_t_graph = function.induced_subgraph(g, visited)
            s_t_graph = nx.Graph(s_t_graph)
            #the induced graph could contain reference edges which need to be removed
            remove = []
            for e in s_t_graph.edges:
                if "%s_0_" % ref_g_id in g.nodes[e[0]]['geneIDs'] \
                and "%s_0_" % ref_g_id in g.nodes[e[1]]['geneIDs']:
                    n1 = mapping.loc[g.nodes[e[0]]["name"]][0]
                    n2 = mapping.loc[g.nodes[e[1]]["name"]][0]
                    if abs(int(n1.split("_")[2]) -
                           int(n2.split("_")[2])) < 100:
                        #print(n1, n2)
                        remove.append(e)
            s_t_graph.remove_edges_from(remove)
            #print some info about that long-range connection
            #print(n)
            #print(nid, sink["sink"])
            #min cut between the two reference nodes
            cut = []
            cut_weight, partitions = nx.algorithms.flow.minimum_cut(
                s_t_graph, nid, sink["sink"])
            for p1_node in partitions[0]:
                for p2_node in partitions[1]:
                    if s_t_graph.has_edge(p1_node, p2_node):
                        cut.append((p1_node, p2_node))
            #cardinality cut TODO make this an option
            #cut = cuts.minimum_edge_cut(s_t_graph, nid, sink["sink"])
            for e in cut:
                print(g.nodes[e[0]]['name'], g.nodes[e[1]]['name'])
                cut_edges.append(e)
            #delete cut edges from the graph
            if len(cut) == 0:
                #something happened as no min cut can be found
                i += 1
                sys.exit(
                    "no min cut could be found; sorry this shouldn't happen")
            g.remove_edges_from(cut)
            sink["sink"] = None
            #there may be more paths from that node -> apply again on the same node
        else:
            #all nodes explored; move on
            i += 1
            sink["sink"] = None
    #write cut edges to disk
    with open(cut_edges_out, "w") as f:
        f.write("shared name\tis_cut_edge\n")
        for e in cut_edges:
            f.write("%s (interacts with) %s\t1\n" %(e[0], e[1]))
            f.write("%s (interacts with) %s\t1\n" %(e[1], e[0]))
    #DEBUG to compress the graph 
    #for n in g.nodes:
    #    gene_ids = [(i.split("_")[0], i) for i in  g.nodes[n]['geneIDs'].split(";")]
    #    gene_ids = list(filter(lambda x: ref_g_id == x[0],gene_ids))
    #    if len(gene_ids) == 1:
    #        g.nodes[n]['geneIDs'] = ""
    #    else:
    #        g.nodes[n]['geneIDs'] = gene_ids[0][1]

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        "enable reference-based layouting through detecting long-range connection between otherwise distant genes in a reference genome"
    )
    parser.add_argument(
        "ref_g_id", help='reference genome id (should be a complete genome)')
    parser.add_argument("graph", help='path to final_graph.gml')
    parser.add_argument("cut_edges_out", help='file for cut edges')
    parser.add_argument(
        "--add_reference_edges",
        action="store_true",
        help=
        'add edges between consecutive genes in the reference genome even if they have been removed by panaroo'
    )
    parser.add_argument("--ignore_high_var",
                        action="store_true",
                        help='ignore highly variable genes')
    args = parser.parse_args()
    layout(**vars(args))
