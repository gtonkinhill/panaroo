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
from scipy import optimize, stats
from numba import jit
from joblib import Parallel, delayed

import pickle

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

            # skip conserved genes
            # if np.sum([int(pa) for pa in line[1:]]) == len(line[1:]):
            #     continue

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
        p = np.logaddexp(v_l - av_l, (a_l - av_l) - (a+v)*t)
    elif (xl==0) and (xn==1):
        p = (a_l - av_l) + log1mexp(-(a+v)*t)
    elif (xl==1) and (xn==0):
        p = (v_l - av_l) + log1mexp(-(a+v)*t)
    else:
        p = np.logaddexp((a_l - av_l), (v_l - av_l) -(a+v)*t)
    return(p)

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

def calc_llk_fmg(params, tree_array, nleaves, presence_absence, isolates):
    
    a = params[0]
    v = params[1]

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

    # calculate llk of conserved genes
    print("Calculating null")
    l0 = np.full(nleaves, -math.inf)
    l1 = np.zeros(nleaves)
    L_conserved = calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v)

    print("L_null: ", L_null)

    print("Calculating llk")
    llk = 0
    for g in tqdm(presence_absence):
        l0 = presence_absence[g][0]
        l1 = presence_absence[g][1]
        llk += calc_llk_gene_numpy(tree_array, nleaves, l0, l1, a, v
            ) - np.log(1 - np.exp(L_null) - np.exp(L_1) - np.exp(L_conserved))

    print("llk: ", llk)

    return llk

@jit(nopython=True) 
def fmg_trans_llk_prob(xl, xn, t, v):
    if (xl==0) and (xn==0):
        p = 0
    elif (xl==0) and (xn==1):
        p = -math.inf
    elif (xl==1) and (xn==0):
        p = log1mexp(-v*t)
    else:
        p = -v*t
    return(p)


def get_origin_nodes(in_tree, node_index, presence_absence, n1=False):
    origin_indices = defaultdict(list)

    if n1:
        for leaf in in_tree.leaf_node_iter():
            mrca = in_tree.mrca(taxon_labels=[leaf.taxon.label])
            origin_nodes = list(mrca.ancestor_iter(inclusive=True))
            for node in origin_nodes:                    
                origin_indices[node_index[leaf.label]].append(node_index[node.label])
    else:
        for gene in presence_absence:
            taxon_labels = [tax for tax in presence_absence[gene] if presence_absence[gene][tax]==1]
            mrca = in_tree.mrca(taxon_labels=taxon_labels)
            origin_nodes = list(mrca.ancestor_iter(inclusive=True))

            for node in origin_nodes:
                origin_indices[gene].append(node_index[node.label])

    return origin_indices

@jit(nopython=True)
def img_llk_all_nodes(in_tree, nleaves, l0, l1, u, v):

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
                        (fmg_trans_llk_prob(xl, xn, in_tree[i][4], v) + 
                        in_tree[int(in_tree[i][0])][2+xn] +
                        fmg_trans_llk_prob(xl, xm, in_tree[i][5], v) + 
                        in_tree[int(in_tree[i][1])][2+xm]))
    
    return(in_tree)

# @jit
def calc_llk_img(params, tree_array, nleaves, origin_nodes, 
    origin_nodes_n1, presence_absence, isolates):
    
    u = params[0]
    v = params[1]

    u_l = np.log(u)
    v_l = np.log(v)

    # calculate expected N_null
    l0 = np.zeros(nleaves)
    l1 = np.full(nleaves, -math.inf)
    tree_array = img_llk_all_nodes(tree_array, nleaves, l0, l1, u, v)
    N_null = -math.inf
    N_tot = -math.inf
    for i in np.arange(nleaves, tree_array.shape[0]):
        if i==(tree_array.shape[0]-1):
            gn_l = u_l - v_l
        else:
            gn_l = u_l - v_l + log1mexp(-v*tree_array[i][6])
        N_null = np.logaddexp(N_null,
            tree_array[i][3] + gn_l)
        N_tot = np.logaddexp(N_tot, gn_l)

    print("N_null: ", N_null)

    # calculate N_1
    N_1 = -math.inf
    l0 = np.zeros(nleaves)
    l1 = np.full(nleaves, -math.inf)
    for i in tqdm(range(nleaves)):    
        l0[i] = -math.inf
        l1[i] = 0
        tree_array = img_llk_all_nodes(tree_array, nleaves, l0, l1, u, v)
        for j in origin_nodes_n1[i]:
            if j==(tree_array.shape[0]-1):
                gn_l = u_l - v_l
            else:
                gn_l = u_l - v_l + log1mexp(-v*tree_array[j][6])
            N_1 = np.logaddexp(N_1,
                tree_array[j][3] + gn_l)
        l0[i] = 0
        l1[i] = -math.inf

    print("N_1: ", N_1)

    # calculate N_total
    for i in np.arange(0, nleaves):
        gn_l = u_l - v_l + log1mexp(-v*tree_array[i][6])
        N_tot = np.logaddexp(N_tot, gn_l)

    print("N_tot: ", N_tot)

    llk = 0
    for g in tqdm(presence_absence):
        N_exp = -math.inf
        l0 = presence_absence[g][0]
        l1 = presence_absence[g][1]
        tree_array = img_llk_all_nodes(tree_array, nleaves, l0, l1, u, v)
        for i in origin_nodes[g]:
            if i==(tree_array.shape[0]-1):
                gn_l = u_l - v_l
            else:
                gn_l = u_l - v_l + log1mexp(-v*tree_array[i][6])
            N_exp = np.logaddexp(N_exp,
                tree_array[i][3] + gn_l)

        llk += N_exp - np.log(np.exp(N_tot)-np.exp(N_exp)-np.exp(N_1))

    print("N_exp: ", N_exp)
    print("llk: ", llk)

    return llk

def get_discrete_gamma_rates(alpha, k):
    points = np.arange(1, 2*k, 2)/(2*k)
    median_rates = stats.gamma.ppf(q=points, a=alpha, scale=1/alpha)
    return(median_rates)


def calc_llk_fmg_with_rate(params, k, tree_array, nleaves, 
    presence_absence, isolates):

    observed_Nall = len(np.unique(presence_absence, axis=0))

    alpha = params[0]
    a_rates = params[1:]
    v_rates = get_discrete_gamma_rates(alpha, k)

    print("alpha:", alpha)
    print("a_rates:", a_rates)
    print("v_rates:", v_rates)

    llk = -math.inf
    for a,v in zip(a_rates, v_rates):
        print(a,v)
        llk = np.logaddexp(llk, calc_llk_fmg((a,v), tree_array, nleaves, 
            presence_absence, isolates) - np.log(k))
    
    return -llk

def calc_llk_img_with_rate(params, k, tree_array, nleaves, origin_nodes, 
    origin_nodes_n1, presence_absence, isolates):

    alpha = params[0]
    u_rates = params[1:]
    v_rates = get_discrete_gamma_rates(alpha, k)

    print("alpha:", alpha)
    print("u_rates:", u_rates)
    print("v_rates:", v_rates)

    llk = -math.inf
    for u,v in zip(u_rates, v_rates):
        print(u,v)
        llk = np.logaddexp(llk, calc_llk_img((u,v), tree_array, nleaves, origin_nodes, 
            origin_nodes_n1, presence_absence, isolates) - np.log(k))
    
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

    # generate array form of tree
    # prepare matrix representation of tree
    nnodes = 0
    for node in tree.leaf_node_iter():
        node.label = nnodes
        nnodes+=1
    for node in tree.postorder_internal_node_iter():
        node.label = nnodes
        nnodes += 1
    tree_array = np.zeros((nnodes, 7))
    leaves = []
    node_index = {}
    for i, node in enumerate(tree.leaf_node_iter()):
        leaves.append(i)
        node_index[node.label] = i
        tree_array[i][0] = -1
        tree_array[i][1] = -1
        tree_array[i][6] = node.edge.length

    nleaves=len(leaves)
    for i, node in enumerate(tree.postorder_internal_node_iter()):
        j=i+nleaves
        node_index[node.label] = j
        children = node.child_nodes()
        
        tree_array[j][0] = node_index[children[0].label]
        tree_array[j][1] = node_index[children[1].label]
        tree_array[j][4] = children[0].edge.length
        tree_array[j][5] = children[1].edge.length
        tree_array[j][6] = node.edge.length


    if args.model=='IMG':

        get origin indices
        print("Obtaining origin nodes...")
        origin_indices = get_origin_nodes(tree, node_index, presence_absence)
        print("Obtaining origin nodes n1...")
        origin_indices_n1 = get_origin_nodes(tree, node_index, presence_absence, n1=True)

        # with open("origin_indices.pkl", 'wb') as outfile:
        #     pickle.dump(origin_indices, outfile)
        # with open("origin_indices_n1.pkl", 'wb') as outfile:
        #     pickle.dump(origin_indices_n1, outfile)

        # with open("origin_indices.pkl", 'rb') as infile:
        #     origin_indices =  pickle.load(infile)
        # with open("origin_indices_n1.pkl", 'rb') as infile:
        #     origin_indices_n1 =  pickle.load(infile)

        alpha_bounds = (1e-2,100)
        u_bounds = (1e-6,1e3)
        k=5
        bounds = [alpha_bounds]
        x0 = [5]
        x0 += [1]*k

        for u in range(k):
            bounds += [u_bounds]

        result = optimize.minimize(calc_llk_img_with_rate, bounds=bounds,
            x0=x0, method='Nelder-Mead',
            args=(k, tree_array, nleaves, origin_indices, origin_indices_n1, 
                presence_absence_llk, isolates))
        print(result)

    else:
        alpha_bounds = (0.01, 1000)
        a_bounds = (1e-6,1000)
        k=5
        bounds = [alpha_bounds]
        x0 = [5]
        x0 += [1]*k

        for u in range(k):
            bounds += [a_bounds]
        result = optimize.minimize(calc_llk_fmg_with_rate, bounds=bounds,
            x0=x0, method='Nelder-Mead',
            args=(k, tree_array, nleaves, presence_absence_llk, isolates))
        print(result)

    return


if __name__ == '__main__':
    main()
