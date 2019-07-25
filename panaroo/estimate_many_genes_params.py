import os
from .__init__ import __version__

import numpy as  np
from dendropy.simulate import treesim
from dendropy.model import reconcile
from dendropy import TaxonNamespace, Tree
import copy
import math
from collections import defaultdict
from tqdm import tqdm
from scipy import optimize
from numba import jit
from joblib import Parallel, delayed

@jit(nopython=True) 
def log1mexp(a):
    if a < np.log(0.5):
        return(np.log1p(-np.exp(a)))
    return(np.log(-np.expm1(a)))

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

            presence_absence[line[0]] = {}
            for iso, pa in zip(isolates, line[1:]):
                presence_absence[line[0]][iso] = int(pa)

    return(isolates, presence_absence)

@jit(nopython=True) 
def trans_llk_prob(xl, xn, t, a, v):
    a_l = np.log(a)
    v_l = np.log(v)
    av_l = np.logaddexp(a_l, v_l)

    if (xl==0) and (xn==0):
        p = np.logaddexp(v_l - av_l, (a_l - av_l) + -(a+v)*t)
    elif (xl==0) and (xn==1):
        p = (a_l - av_l) + log1mexp(-(a+v)*t)
    elif (xl==1) and (xn==0):
        p = (v_l - av_l) + log1mexp(-(a+v)*t)
    else:
        p = (a_l - av_l) + (v_l - av_l) -(a+v)*t

    return(p)

# def calc_llk_gene(in_tree, presence_absence, a, v):
    
#     # set lead nodes
#     for node in in_tree.leaf_node_iter():
#         node.llk = [0.0, 0.0]
#         node.llk[0] = 0.0 if presence_absence[node.taxon.label]==0 else -math.inf
#         node.llk[1] = 0.0 if presence_absence[node.taxon.label]==1 else -math.inf

#     # now iterate bottom up to fill in the values for internal nodes
#     for node in in_tree.postorder_internal_node_iter():
#         node.llk = [-math.inf, -math.inf]
#         for xl in [0,1]:
#             for xn in [0,1]:
#                 for xm in [0,1]:
#                     children = node.child_nodes()                    
#                     node.llk[xl] = np.logaddexp(node.llk[xl],
#                         (trans_llk_prob(xl, xn, children[0].edge.length, a, v) + 
#                         children[0].llk[xn] +
#                         trans_llk_prob(xl, xm, children[1].edge.length, a, v) + 
#                         children[1].llk[xm]))

#     llk = np.logaddexp(in_tree.seed_node.llk[1]+np.log(a)-np.log(a+v),
#         in_tree.seed_node.llk[0]+np.log(v)-np.log(a+v))

#     return(llk)

@jit(nopython=True) 
def calc_llk_gene_numpy(in_tree, nleaves, l0, l1, a, v):

    # set lead nodes
    in_tree[0:nleaves,2] = l0
    in_tree[0:nleaves,3] = l1

    for i in np.arange(nleaves, in_tree.shape[0]):
        in_tree[i][2] = -math.inf
        in_tree[i][3] = -math.inf
        for xl in [0,1]:
            for xn in [0,1]:
                for xm in [0,1]:
                    in_tree[i][2+xl] = np.logaddexp(in_tree[i][2+xl],
                        (trans_llk_prob(xl, xn, in_tree[i][4], a, v) + 
                        in_tree[int(in_tree[i][0])][2+xn] +
                        trans_llk_prob(xl, xm, in_tree[i][5], a, v) + 
                        in_tree[int(in_tree[i][1])][2+xm]))

    llk = np.logaddexp(in_tree[in_tree.shape[0]-1][3]+np.log(a)-np.log(a+v),
        in_tree[in_tree.shape[0]-1][2]+np.log(v)-np.log(a+v))

    return(llk)

def calc_llk(params, in_tree, presence_absence, isolates):
    
    a = params[0]
    v = params[1]
    print("a:", a, "  v:", v)
    # prepare matrix representation of tree
    nnodes = 0
    for node in in_tree.leaf_node_iter():
        node.label = nnodes
        nnodes+=1
    for node in in_tree.postorder_internal_node_iter():
        node.label = nnodes
        nnodes += 1
    tree_array = np.zeros((nnodes, 6))
    leaves = []
    node_index = {}
    for i, node in enumerate(in_tree.leaf_node_iter()):
        leaves.append(i)
        node_index[node.label] = i
        tree_array[i][0] = -1
        tree_array[i][1] = -1

    nleaves=len(leaves)
    for i, node in enumerate(in_tree.postorder_internal_node_iter()):
        j=i+nleaves
        node_index[node.label] = j
        children = node.child_nodes()
        tree_array[j][0] = node_index[children[0].label]
        tree_array[j][1] = node_index[children[1].label]
        tree_array[j][4] = children[0].edge.length
        tree_array[j][5] = children[1].edge.length

    # calculate llk of null pattern
    print("Calculating null")
    l0 = np.zeros(nleaves)
    l1 = np.full(nleaves, -math.inf)
    L_null = calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v)

    print("L_null: ", L_null)

    # calculate llk of combined L1 pattern
    print("Calculating L1")
    L_1 = -math.inf
    l0 = np.zeros(nleaves)
    l1 = np.full(nleaves, -math.inf)
    for i in tqdm(range(nleaves)):    
        l0[i] = -math.inf
        l1[i] = 0
        L_1 = np.logaddexp(L_1,
            calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v))
        l0[i] = 0
        l1[i] = -math.inf

    print("L_1: ", L_1)

    print("Calculating llk")
    llk = -math.inf
    for g in tqdm(presence_absence):
        l0 = presence_absence[g][0]
        l1 = presence_absence[g][1]
        llk = np.logaddexp(llk, calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v) - 
            np.log(1 - np.exp(L_null) - np.exp(L_1)))

    print("llk: ", llk)

    return -llk


def get_options():
    import argparse

    description = 'Estimate model parameters for either the Infinitely or Finitely Many Genes Model'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo_mg_est')

    parser.add_argument(
        "--model",
        dest="model",
        help=
        "Specify a model to estimate. One of 'IMG' or 'FMG'",
        type=str,
        choices={'IMG', 'FMG'},
        default="IMG")

    parser.add_argument(
        "--tree",
        dest="tree",
        help="A dated phylogeny.",
        type=str,
        required=True)

    parser.add_argument(
        "--pa",
        dest="presence_absence",
        help="A presence/absence produced by Panaroo.",
        type=str,
        required=True)

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

    # gather gene prescence abscence for convinience in llk calculation
    presence_absence_llk = {}
    for gene in presence_absence:
        presence_absence_llk[gene] = [np.array([0.0 if 
            presence_absence[gene][node.taxon.label]==0 else 
            -math.inf for node in tree.leaf_node_iter()]),
            np.array([0.0 if presence_absence[gene][node.taxon.label]==1 
                else -math.inf for node in tree.leaf_node_iter()])]

    calc_llk((1e-6,5), tree, presence_absence_llk, isolates)

    # find some decent paramters
    # bounds = ((1e-6,5), (1e-6,5))
    # result = optimize.shgo(calc_llk, bounds=bounds,
    #     args=(tree, presence_absence_llk, isolates))
    # print(result)

    return


if __name__ == '__main__':
    main()
