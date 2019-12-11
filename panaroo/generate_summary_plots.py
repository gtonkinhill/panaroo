import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt


def plot_isolates_per_gene_hist(G):

    plt.style.use('ggplot')

    node_gene_nos = []

    for node in G.nodes():
        node_gene_nos.append(len(G.nodes[node]["members"]))

    fig = plt.figure()
    axis1 = fig.add_subplot(111)
    axis1.hist(node_gene_nos, bins=20)
    axis1.set_xlabel("Number of genes")
    axis1.set_ylabel("Number of isolates")
    fig.savefig("isolates_per_gene_hist.png")

    return
