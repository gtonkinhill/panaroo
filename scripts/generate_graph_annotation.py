import networkx as nx

def generate_annotation(annotation, graph, node_out, edge_out):
    annot_dicts = {} 
    #resolve hzi isolate name
    #read in graph and get isolate names
    g = nx.read_gml(graph)
    import re
    #id2isol = dict([(re.sub(r'_spades3.10.0_careful$', '', i), str(j)) for i, j in zip(g.graph["isolateNames"], range(len(g.graph["isolateNames"])))]) #resolve hzi isolate name clash
    id2isol = dict([(re.sub(r'\.velvet', '', i), str(j)) for i, j in zip(g.graph["isolateNames"], range(len(g.graph["isolateNames"])))]) #resolve hzi isolate name clash
    #id2isol = dict([(i, j) for i, j in zip(g.graph["isolateNames"], range(len(g.graph["isolateNames"])))])
    #parse user annotation
    with open(annotation, 'r') as annot:
        #header
        header = [i for i in annot.readline().strip().split("\t")]
        header = header[1:]
        field2annot = dict(zip(range(len(header)), header))
        for i in header:
            annot_dicts[i] = {}
        for l in annot:
            fields = l.strip().split("\t")
            #check if isolate is not in the graph and ignore if this is the case
            if fields[0] not in id2isol:
                continue
            isol = id2isol[fields[0]]
            fields = fields[1:]
            for f in range(len(field2annot)):
                annot_dicts[field2annot[f]][isol] = fields[f]
    annot_node_dicts = {}
    #set up node dictionary for each annotation column
    for annot in annot_dicts:
        annot_node_dicts[annot] = {} 
    for n in g.nodes:
        name = g.nodes[n]['name']
        for m in g.nodes[n]["members"]:
            if m not in annot_dicts[annot]:
                continue
            for annot in annot_node_dicts:
                if name not in annot_node_dicts[annot]:
                    annot_node_dicts[annot][name] = set([annot_dicts[annot][m]])
                else:
                    annot_node_dicts[annot][name].add(annot_dicts[annot][m])
    #get all nodes
    nodes = [g.nodes[n]['name'] for n in g.nodes]
    with open(node_out, 'w') as no:
        no.write("shared name\t")
        no.write("\t".join(list(annot_node_dicts.keys())) + "\n")
        for n in nodes:
            no.write(n)
            for annot in annot_node_dicts:
                no.write("\t")
                if n in annot_node_dicts[annot]:
                    no.write(";".join(annot_node_dicts[annot][n]))
            no.write("\n") 
    annot_edge_dicts = {}
    for annot in annot_dicts:
        annot_edge_dicts[annot] = {} 
    for e in g.edges:
        s_t = e[0] + "_" + e[1] 
        for m in g.edges[e]["members"]:
            if m not in annot_dicts[annot]:
                continue
            for annot in annot_edge_dicts:
                if s_t not in annot_edge_dicts[annot]:
                    annot_edge_dicts[annot][s_t] = set([annot_dicts[annot][m]])
                else:
                    annot_edge_dicts[annot][s_t].add(annot_dicts[annot][m])

    edges = list(annot_edge_dicts[list(annot_edge_dicts.keys())[0]].keys())
    with open(edge_out, 'w') as eo:
        eo.write("shared name" + "\t")
        eo.write("\t".join(list(annot_edge_dicts.keys())) + "\n")
        for e in edges:
            eo.write(" (interacts with) ".join(e.split("_")))
            for annot in annot_edge_dicts:
                eo.write("\t")
                if e in annot_edge_dicts[annot]:
                    eo.write(";".join(annot_edge_dicts[annot][e]))
            eo.write("\n") 

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        "generate graph annotation from user provided isolate annotation e.g. isolate sequence type or AMR information"
    )
    parser.add_argument(
        "annotation", help='table with one annotation type per column')
    parser.add_argument(
        "graph", help='final graph panaroo output')
    parser.add_argument(
        "edge_out", help='edge output mapping file')
    parser.add_argument(
        "node_out", help='node output mapping file')
    args = parser.parse_args()
    generate_annotation(**vars(args))
