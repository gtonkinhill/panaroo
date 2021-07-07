import os
from .__init__ import __version__

import numpy as np
from dendropy.simulate import treesim
from dendropy.model import reconcile
from dendropy import TaxonNamespace, Tree
import copy
import math
from collections import defaultdict
from scipy import optimize, stats
from numba import jit
from joblib import Parallel, delayed
import random
from scipy.stats.distributions import chi2


@jit(nopython=True)
def log1mexp(a):
    if a < np.log(0.5):
        return (np.log1p(-np.exp(-a)))
    return (np.log(-np.expm1(-a)))


@jit(nopython=True)
def log_subtract(x, y):
    if x <= y:
        raise RuntimeError("error!! computing the log of a negative number")
    if y == -np.inf:
        return (x)
    return (x + np.log1p(-np.exp(y - x)))


def load_pa(presence_absence_file):

    with open(presence_absence_file, 'r') as infile:
        header = next(infile)
        presence_absence = {}
        isolates = header.strip().split()[1:]
        isolates = [iso.replace("_", " ") for iso in isolates]

        for line in infile:
            line = line.strip().split()

            # skip singleton genes
            if np.sum([int(pa) for pa in line[1:]]) <= 1:
                continue

            # skip conserved genes
            # if np.sum([int(pa) for pa in line[1:]]) == len(line[1:]):
            #     continue

            presence_absence[line[0]] = {}
            for iso, pa in zip(isolates, line[1:]):
                presence_absence[line[0]][iso] = int(pa)

    return (isolates, presence_absence)


@jit(nopython=True)
def trans_llk_prob(xl, xn, t, a, v):
    a_l = np.log(a)
    v_l = np.log(v)
    av_l = np.logaddexp(a_l, v_l)
    if (xl == 0) and (xn == 0):
        p = np.logaddexp(v_l - av_l, (a_l - av_l) - (a + v) * t)
    elif (xl == 0) and (xn == 1):
        p = (a_l - av_l) + log1mexp((a + v) * t)
    elif (xl == 1) and (xn == 0):
        p = (v_l - av_l) + log1mexp((a + v) * t)
    else:
        p = np.logaddexp((a_l - av_l), (v_l - av_l) - (a + v) * t)
    return (p)


@jit(nopython=True)
def calc_llk_gene_numpy(in_tree, nleaves, l0, l1, a, v):

    # set lead nodes
    in_tree[0:nleaves, 2] = l0
    in_tree[0:nleaves, 3] = l1

    for i in np.arange(nleaves, in_tree.shape[0]):
        in_tree[i][2] = -math.inf
        in_tree[i][3] = -math.inf
        for xl in [0, 1]:
            for xn in [0, 1]:
                for xm in [0, 1]:
                    in_tree[i][2 + xl] = np.logaddexp(
                        in_tree[i][2 + xl],
                        (trans_llk_prob(xl, xn, in_tree[i][4], a, v) +
                         in_tree[int(in_tree[i][0])][2 + xn] +
                         trans_llk_prob(xl, xm, in_tree[i][5], a, v) +
                         in_tree[int(in_tree[i][1])][2 + xm]))

    llk = np.logaddexp(
        in_tree[in_tree.shape[0] - 1][3] + np.log(a) - np.log(a + v),
        in_tree[in_tree.shape[0] - 1][2] + np.log(v) - np.log(a + v))

    return (llk)


def calc_llk_fmg(params, tree_array, nleaves, presence_absence, isolates,
                 verbose):

    a = params[0]
    v = params[1]

    if verbose: print("a:", a, " v:", v)

    if (a < 0) or (v < 0):
        return (-math.inf)

    # calculate llk of null pattern
    if verbose: print("Calculating null")
    l0 = np.zeros(nleaves)
    l1 = np.full(nleaves, -math.inf)
    L_null = calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v)

    if verbose: print("L_null: ", L_null)

    # calculate llk of combined L1 pattern
    if verbose: print("Calculating L1")
    L_1 = -math.inf
    l0 = np.zeros(nleaves)
    l1 = np.full(nleaves, -math.inf)
    for i in range(nleaves):
        l0[i] = -math.inf
        l1[i] = 0
        L_1 = np.logaddexp(
            L_1, calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v))
        l0[i] = 0
        l1[i] = -math.inf
    if verbose: print("L_1: ", L_1)

    # calculate llk of conserved genes
    # if verbose: print("Calculating null")
    # l0 = np.full(nleaves, -math.inf)
    # l1 = np.zeros(nleaves)
    # L_conserved = calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v)

    # if verbose: print("L_conserved: ", L_conserved)

    if verbose: print("Calculating llk")
    llk = 0
    for g in presence_absence:
        l0 = presence_absence[g][0]
        l1 = presence_absence[g][1]
        llk += calc_llk_gene_numpy(
            tree_array, nleaves, l0, l1, a,
            v) - np.log(1 - np.exp(L_null) - np.exp(L_1))
        # print(L_null, L_1)
        # print(1 - np.exp(L_null) - np.exp(L_1))
        # print("llka: ", calc_llk_gene_numpy(
        #     tree_array, nleaves, l0, l1, a,
        #     v))
        # print("llkb: ", llk)

    if verbose: print("llk: ", llk)

    return -llk


def get_discrete_gamma_rates(alpha, k):
    points = np.arange(1, 2 * k, 2) / (2 * k)
    median_rates = stats.gamma.ppf(q=points, a=alpha, scale=1 / alpha)
    return (median_rates)


def get_options():
    import argparse

    description = 'Estimate model parameters for either the Finitely Many Genes Model'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo_fmg_est')

    parser.add_argument("--tree",
                        dest="tree",
                        help="A dated phylogeny.",
                        type=str,
                        required=True)

    parser.add_argument("--pa",
                        dest="presence_absence",
                        help="A presence/absence produced by Panaroo.",
                        type=str,
                        required=True)

    parser.add_argument("-o",
                        "--outfile",
                        dest="outputfile",
                        help="Name of outputfile.",
                        type=str,
                        required=True)

    parser.add_argument(
        "--nboot",
        dest="nboot",
        help="The number of sub-sampling bootstrap iterations to perform.",
        type=int,
        default=0)

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    parser.add_argument("--verbose",
                        dest="verbose",
                        help="print additional output",
                        action='store_true',
                        default=False)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)


def optimise_model(model,
                   presence_absence_llk,
                   bounds,
                   x0,
                   tree_array,
                   nleaves,
                   isolates,
                   origin_nodes=None,
                   origin_nodes_n1=None):
    all_genes = list(presence_absence_llk.keys())
    boot_pa = {}
    new_genes = random.choices(all_genes, k=len(all_genes))
    boot_gene_count = np.zeros(len(isolates))
    for j, g in enumerate(new_genes):
        boot_pa[j] = presence_absence_llk[g]
        boot_gene_count += presence_absence_llk[g][1] == 0

    boot_result = optimize.minimize(calc_llk_fmg,
                                    bounds=bounds,
                                    x0=x0,
                                    method='L-BFGS-B',
                                    args=(tree_array, nleaves, boot_pa,
                                          isolates, False))

    return ((boot_result.x[0], boot_result.x[1], np.mean(boot_gene_count)))


def main():
    args = get_options()

    # load input data
    isolates, presence_absence = load_pa(args.presence_absence)
    tree = Tree.get(path=args.tree, schema="newick")

    # check tree nodes match isolates
    leaf_taxa = set()
    for node in tree.leaf_node_iter():
        node.taxon.label = node.taxon.label.strip()
        node.taxon.label = node.taxon.label.strip("'")
        node.taxon.label = node.taxon.label.strip('"')
        leaf_taxa.add(node.taxon.label)
        if node.taxon.label not in isolates:
            raise RuntimeError("Mismatch between tree and p/a matrix!")
    for iso in isolates:
        if iso not in leaf_taxa:
            raise RuntimeError("Mismatch between tree and p/a matrix!")

    # gather gene prescence abscence for convenience in llk calculation
    presence_absence_llk = {}
    for gene in presence_absence:
        presence_absence_llk[gene] = [
            np.array([
                0.0
                if presence_absence[gene][node.taxon.label] == 0 else -math.inf
                for node in tree.leaf_node_iter()
            ]),
            np.array([
                0.0
                if presence_absence[gene][node.taxon.label] == 1 else -math.inf
                for node in tree.leaf_node_iter()
            ])
        ]

    # generate array form of tree
    # prepare matrix representation of tree
    nnodes = 0
    for node in tree.leaf_node_iter():
        node.label = nnodes
        nnodes += 1
    for node in tree.postorder_internal_node_iter():
        node.label = nnodes
        nnodes += 1
    tree_array = np.zeros((nnodes, 7))
    leaves = []
    node_index = {}
    for i, node in enumerate(tree.leaf_node_iter()):
        if node.edge.length <= 0: 
            print('edge length: ', node.edge.length)
            raise RuntimeError('Tree edge length must > 0!')
        leaves.append(i)
        node_index[node.label] = i
        tree_array[i][0] = -1
        tree_array[i][1] = -1
        tree_array[i][6] = node.edge.length

    nleaves = len(leaves)
    for i, node in enumerate(tree.postorder_internal_node_iter()):
        if node.edge.length <= 0: 
            print('edge length: ', node.edge.length)
            raise RuntimeError('Tree edge length must > 0!')
        j = i + nleaves
        node_index[node.label] = j
        children = node.child_nodes()

        tree_array[j][0] = node_index[children[0].label]
        tree_array[j][1] = node_index[children[1].label]
        tree_array[j][4] = children[0].edge.length
        tree_array[j][5] = children[1].edge.length
        tree_array[j][6] = node.edge.length

    for i in range(tree_array.shape[0]):
        print(tree_array[i][6])

    outfile = open(args.outputfile, 'w')

    a_bounds = (1e-7, 1e3)
    v_bounds = (1e-7, 1e3)
    bounds = [a_bounds, v_bounds]
    x0 = [0.001, 0.001]

    result = optimize.minimize(calc_llk_fmg,
                               bounds=bounds,
                               x0=x0,
                               method='L-BFGS-B',
                               args=(tree_array, nleaves, presence_absence_llk,
                                     isolates, args.verbose))

    gene_count = np.zeros(len(isolates))
    for g in presence_absence_llk:
        gene_count += presence_absence_llk[g][1] == 0
    gene_count = np.mean(gene_count)

    if args.nboot > 0:
        # all_genes = list(presence_absence_llk.keys())
        all_boots = np.zeros((args.nboot, 3))

        boot_results = Parallel(n_jobs=args.n_cpu)(
            delayed(optimise_model)("FMG", presence_absence_llk, bounds, x0,
                                    tree_array, nleaves, isolates)
            for i in range(args.nboot))
        for i, boot in enumerate(boot_results):
            all_boots[i, 0] = boot[0]
            all_boots[i, 1] = boot[1]
            all_boots[i, 2] = boot[2]

        a_q = np.quantile(all_boots[:, 0], np.array([0.025, 0.975]))
        v_q = np.quantile(all_boots[:, 1], np.array([0.025, 0.975]))
        G_q = np.quantile(all_boots[:, 2], np.array([0.025, 0.975]))
    else:
        a_q = [math.nan, math.nan]
        v_q = [math.nan, math.nan]
        G_q = [math.nan, math.nan]

    outfile.write("Model: Finitely Many Genes (FMG)\n")
    outfile.write("Llk: " + str(-result.fun) + "\n")
    outfile.write(f"BIC: {2*np.log(len(isolates)) + 2*result.fun}\n")
    outfile.write("Parameter,2.5%CI,97.5%CI\n")
    outfile.write(f"a,{result.x[0]},{a_q[0]},{a_q[1]}\n")
    outfile.write(f"v,{result.x[1]},{v_q[0]},{v_q[1]}\n")
    outfile.write(f"G,{gene_count},{G_q[0]},{G_q[1]}\n")
    outfile.write(
        f"M,{gene_count*(result.x[0]+result.x[1])/result.x[0]},{G_q[0]*(a_q[0]+v_q[0])/a_q[0]},{G_q[1]*(a_q[1]+v_q[1])/a_q[1]}\n"
    )

    outfile.close()

    return


if __name__ == '__main__':
    main()
