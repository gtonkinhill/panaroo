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

def load_pa(presence_absence_file):

    with open(presence_absence_file, 'r') as infile:
        header = next(infile)
        presence_absence = {}
        isolates = header.strip().split()[1:]
        isolates = [iso.replace("_", " ") for iso in isolates]

        for line in infile:
            line = line.strip().split()
            presence_absence[line[0]] = {}
            for iso, pa in zip(isolates, line[1:]):
                presence_absence[line[0]][iso] = int(pa)

    return(isolates, presence_absence)

def trans_llk_prob(xl, xn, t, a, v):
    
    if (xl==0) and (xn==0):
        p = v/(a+v)+(a/(a+v))*np.exp(-(a+v)*t)
    elif (xl==0) and (xn==1):
        p = (a/(a+v))*(1-np.exp(-(a+v)*t))
    elif (xl==1) and (xn==0):
        p = (v/(a+v))*(1-np.exp(-(a+v)*t))
    else:
        p = a/(a+v)+(v/(a+v))*np.exp(-(a+v)*t)

    return(np.log(p))


def calc_llk_gene(in_tree, presence_absence, a, v):
    # set lead nodes
    for node in in_tree.leaf_node_iter():
        node.llk = [0.0, 0.0]
        node.llk[0] = 0.0 if presence_absence[node.taxon.label]==0 else -math.inf
        node.llk[1] = 0.0 if presence_absence[node.taxon.label]==1 else -math.inf

    # now iterate bottom up to fill in the values for internal nodes
    for node in in_tree.postorder_internal_node_iter():
        node.llk = [0.0, 0.0]
        for xl in [0,1]:
            for xn in [0,1]:
                for xm in [0,1]:
                    children = node.child_nodes()                    
                    node.llk[xl] = np.logaddexp(node.llk[xl],
                        (trans_llk_prob(xl, xn, children[0].edge.length, a, v) + 
                        children[0].llk[xn] +
                        trans_llk_prob(xl, xm, children[1].edge.length, a, v) + 
                        children[1].llk[xm]))

    llk = np.logaddexp(in_tree.seed_node.llk[1]+np.log(a)-np.log(a+v),
        in_tree.seed_node.llk[0]+np.log(v)-np.log(a+v))

    return(llk)


def calc_llk(in_tree, presence_absence, isolates, a, v):
    
    # calculate llk of null pattern
    print("Calculating null")
    null_presence_absence = {}
    for tax in isolates:
        null_presence_absence[tax] = 0
    L_null = calc_llk_gene(in_tree, null_presence_absence, a, v)

    # calculate llk of combined L1 pattern
    print("Calculating L1")
    L_1 = 0.0
    for tax in tqdm(isolates):
        null_presence_absence[tax] = 1
        L_1 = np.logaddexp(L_1,
            calc_llk_gene(in_tree, null_presence_absence, a, v))
        null_presence_absence[tax] = 0

    print("Calculating llk")
    llk = 0
    for g in tqdm(presence_absence):
        llk = np.logaddexp(llk, calc_llk_gene(in_tree, presence_absence[g], a, v) - 
            np.log(1 - np.exp(L_null) - np.exp(L_1)))

    return llk


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

    calc_llk(tree, presence_absence, isolates, 1, 1)

    return


if __name__ == '__main__':
    main()
