import numpy as np
import os
import argparse
from .isvalid import *
from .__init__ import __version__
import dendropy as dpy
from collections import defaultdict, Counter
from tqdm import tqdm

def read_presence_absence(filename):
    # load matrix into binary array
    matrix_txt = np.loadtxt(filename, dtype=str, delimiter="\t", comments=None)
    sample_names = matrix_txt[0, 1:]
    gene_names = matrix_txt[1:, 0]
    pa_matrix = matrix_txt[1:, 1:] == "1"

    # remove conserved genes
    keep = np.sum(pa_matrix, axis=1) != len(sample_names)
    pa_matrix = pa_matrix[keep, :]
    gene_names = gene_names[keep]

    return (pa_matrix, gene_names, sample_names)


def get_weights_phylogeny(tree_file, sample_names):

    # read in tree
    tree = dpy.Tree.get(path=tree_file,
                        schema="newick",
                        preserve_underscores=True)

    # check sample names match up
    tip_labels = [t.label for t in tree.taxon_namespace]

    if len(set(tip_labels).intersection(set(sample_names))) != len(tip_labels):
        raise ValueError(
            "Sample name mismatch between tree file and presence/absence matrix"
        )

    # iterate top down to weight by edge length.
    for node in tree.preorder_node_iter():
        if node.parent_node is not None:
            node.total_weight = node.parent_node.total_weight + float(
                node.edge.length) / len(node.leaf_nodes())
        else:
            node.total_weight = 0

    # iterate through leaves getting weights
    weight_dict = {}
    for node in tree.leaf_node_iter():
        weight_dict[node.taxon.label] = node.total_weight

    weights = np.zeros(len(sample_names))
    for i, name in enumerate(sample_names):
        weights[i] = weight_dict[name]

    # normalise (not really necessary)
    weights = weights / np.sum(weights)

    return weights


def get_weights_cluster_csv(cluster_file, sample_names):

    weight_dict = {}
    cluster_count = Counter()
    with open(cluster_file, 'r') as infile:
        for line in infile:
            line = line.strip().split(",")
            weight_dict[line[0]] = line[1]
            cluster_count[line[1]] += 1

    # now check sample names match up
    weights = np.zeros(len(sample_names))
    for i, name in enumerate(sample_names):
        if name not in weight_dict:
            raise ValueError(
                "Sample name mismatch between weights and presence/absence matrix"
            )
        weights[i] = 1.0/cluster_count[weight_dict[name]]

    return weights


def spydrpick(pa_matrix, weights=None, keep_quantile=0.9, chunk_size=100):

    ngenes = pa_matrix.shape[0]
    nsamples = pa_matrix.shape[1]

    # set weights to equal if none are specified.
    if weights is None:
        weights = np.ones(nsamples)

    if len(weights) != nsamples:
        raise ValueError("Weight/matrix dimension mismatch!")

    n_eff = np.sum(weights)

    # determine threshold for storing results
    sample_rows = np.random.randint(0, ngenes, size=min(100, ngenes))
    mi_11 = (np.inner(weights * pa_matrix[sample_rows, :], pa_matrix) +
             0.5) / (n_eff + 2.0)
    mi_10 = (np.inner(weights * pa_matrix[sample_rows, :], 1 - pa_matrix) +
             0.5) / (n_eff + 2.0)
    mi_01 = (np.inner(weights * (1 - pa_matrix[sample_rows, :]), pa_matrix) +
             0.5) / (n_eff + 2.0)
    mi_00 = (np.inner(weights *
                      (1 - pa_matrix[sample_rows, :]), 1 - pa_matrix) +
             0.5) / (n_eff + 2.0)

    mi_1y = mi_11 + mi_10
    mi_0y = mi_01 + mi_00
    mi_x1 = mi_11 + mi_01
    mi_x0 = mi_10 + mi_00

    mi = mi_00 * (np.log(mi_00) - np.log(mi_0y) - np.log(mi_x0))
    mi += mi_01 * (np.log(mi_01) - np.log(mi_0y) - np.log(mi_x1))
    mi += mi_10 * (np.log(mi_10) - np.log(mi_1y) - np.log(mi_x0))
    mi += mi_11 * (np.log(mi_11) - np.log(mi_1y) - np.log(mi_x1))

    threshold = np.quantile(mi, keep_quantile)

    # operate in chunks to speed things up
    hitsA = []
    hitsB = []
    mis = []

    for i in range(0, ngenes, chunk_size):
        # calculate co-occurence matrix using matrix algebra
        mi_11 = (np.inner(weights * pa_matrix[i:(i + 100), :], pa_matrix) +
                 0.5) / (n_eff + 2.0)
        mi_10 = (np.inner(weights * pa_matrix[i:(i + 100), :], 1 - pa_matrix) +
                 0.5) / (n_eff + 2.0)
        mi_01 = (np.inner(weights *
                          (1 - pa_matrix[i:(i + 100), :]), pa_matrix) +
                 0.5) / (n_eff + 2.0)
        mi_00 = (np.inner(weights *
                          (1 - pa_matrix[i:(i + 100), :]), 1 - pa_matrix) +
                 0.5) / (n_eff + 2.0)

        mi_1y = mi_11 + mi_10
        mi_0y = mi_01 + mi_00
        mi_x1 = mi_11 + mi_01
        mi_x0 = mi_10 + mi_00

        mi = mi_00 * (np.log(mi_00) - np.log(mi_0y) - np.log(mi_x0))
        mi += mi_01 * (np.log(mi_01) - np.log(mi_0y) - np.log(mi_x1))
        mi += mi_10 * (np.log(mi_10) - np.log(mi_1y) - np.log(mi_x0))
        mi += mi_11 * (np.log(mi_11) - np.log(mi_1y) - np.log(mi_x1))

        np.fill_diagonal(mi, threshold - 1)

        hitA, hitB = np.where(mi > threshold)
        mi = mi[hitA, hitB]
        hitA += i
        keep = hitA!=hitB

        # append current chunk output
        hitsA.append(hitA[keep])
        hitsB.append(hitB[keep])
        mis.append(mi[keep])

    hits = set()
    for a,b,m in zip(np.concatenate(hitsA, axis=0), 
                    np.concatenate(hitsB, axis=0),
                    np.concatenate(mis, axis=0)):
        if a<=b:
            hits.add((a,b,m))
        else:
            hits.add((b,a,m))

    return hits


def tukey_outlier(hitsA, hitsB, mis):

    ids = np.unique(hitsA)
    max_hit_mis = np.zeros(len(ids))

    for i in range(len(ids)):
        max_hit_mis[i] = np.max(mis[hitsA == ids[i]])

    Q1, Q3 = np.quantile(max_hit_mis, [0.25, 0.75])

    outliers = np.zeros(len(mis))
    outliers[mis > (Q3 + 1.5 * (Q3 - Q1))] = 1
    outliers[mis > (Q3 + 3 * (Q3 - Q1))] = 2

    return outliers


def aracne(hitsA, hitsB, mis):
    # build structures to speed up filtering
    allg = set(hitsA) | set(hitsB)

    d = {}
    for a,b,m in zip(hitsA, hitsB, mis):
        d[(a,b)] = m
        d[(b,a)] = m

    # aracne algorithm
    bad_pairs = set()
    nhitsA = []
    nhitsB = []
    nmis = []
    for a,b,m in tqdm(zip(hitsA, hitsB, mis)):
        bad = False
        for c in allg:
            if c in [a,b]: continue
            if (a, c) not in d: continue
            if (b, c) not in d: continue
            if (m < d[(a, c)]) and (m < d[(b, c)]):
                bad = True
        if not bad:
            nhitsA.append(a)
            nhitsB.append(b)
            nhitsA.append(m)

    return (nhitsA, nhitsB, nmis)


def get_options():
    import argparse

    description = 'Identify associations between genes using the Spydrpick algorithm'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo-spydrpick')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-i",
        "--input",
        dest="pa_file",
        required=True,
        help="gene presence absence file (.Rtab) generated by Panaroo")
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=lambda x: is_valid_folder(parser, x))

    io_opts.add_argument(
        "--tree",
        dest="tree_file",
        default=None,
        help=("phylogeny in newick format for weigting samples to" +
              " control for population structure"))

    io_opts.add_argument(
        "--clusters",
        dest="cluster_file",
        default=None,
        help=("sample clusters for weigting to control for " +
              "population structure. format: 'sample_name,cluster_id'"))

    # Other options
    parser.add_argument(
        "--quantile",
        dest="quantile",
        default=0.9,
        type=float,
        help=
        "the quantile used to determine a threshold for keeping MI values (default=0.9)."
    )

    parser.add_argument(
        "--aracne",
        dest="aracne",
        default=False,
        action='store_true',
        help=
        "turn on the Aracne stage of the algorithm (see the Spydrpick paper for details)."
    )

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

    if args.tree_file is not None:
        if args.cluster_file is not None:
            raise ValueError(
                "Only one of 'tree' or 'clusters' may be used for weighting!")

    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    pa_matrix, gene_names, sample_names = read_presence_absence(args.pa_file)

    # set weights
    if args.tree_file is not None:
        weights = get_weights_phylogeny(args.tree_file, sample_names)
    elif args.cluster_file is not None:
        weights = get_weights_cluster_csv(args.cluster_file, sample_names)
    else:
        weights = None
    
    hitsA, hitsB, mis = zip(*spydrpick(pa_matrix,
                                  weights=weights,
                                  keep_quantile=args.quantile,
                                  chunk_size=100))

    outliers = tukey_outlier(hitsA, hitsB, mis)

    if args.aracne:
        hitsA, hitsB, mis = aracne(hitsA, hitsB, mis)

    with open(args.output_dir + "gene_pa_spydrpick.csv", 'w') as outfile:
        outfile.write("GeneA,GeneB,MI,outlier\n")
        for i in np.argsort(-mis):
            outfile.write(",".join([
                gene_names[hitsA[i]], gene_names[hitsB[i]],
                str(mis[i]),
                str(outliers[i])
            ]) + "\n")

    return


if __name__ == '__main__':
    main()
