import os
from .__init__ import __version__

import numpy as np
from dendropy.simulate import treesim
from dendropy.model import reconcile
from dendropy import TaxonNamespace, Tree
import copy
import math
from collections import defaultdict, Counter
from tqdm import tqdm
from scipy import optimize, stats
from .isvalid import *


def log1mexp(a):
    if a < np.log(0.5):
        return (np.log1p(-np.exp(a)))
    return (np.log(-np.expm1(a)))


def load_pa(presence_absence_file):

    with open(presence_absence_file, 'r') as infile:
        header = next(infile)
        isolates = header.strip().split()[1:]
        isolates = [iso.replace("_", " ") for iso in isolates]
        genes = []
        for line in infile:
            line = line.strip().split()
            genes.append(line[0])

    pa_matrix = np.loadtxt(presence_absence_file,
                           delimiter="\t",
                           skiprows=1,
                           usecols=range(1,
                                         len(isolates) + 1),
                           dtype=int)

    return (isolates, genes, pa_matrix)


def f_getspectrum(pa_matrix):
    gene_sizes = np.bincount(np.sum(pa_matrix, 1),
                             minlength=pa_matrix.shape[1] + 1)
    Gk = gene_sizes[1:]  #ingore count of zeros
    return (Gk)


def f_meanpancore(Gk):

    ngenomes = len(Gk)

    gk_pan = np.zeros(ngenomes)  #mean pangenome curve from G(k)
    gk_core = np.zeros(ngenomes)  #mean core curve from G(k)

    #a gene is present in k genomes out of ng
    #the probability that the gene is absent in n genomes out of ng is...
    Pabs = np.zeros((ngenomes, ngenomes))
    #the probability that the gene is present in all n genomes out of ng is...
    Pall = np.zeros((ngenomes, ngenomes))

    #calculate Pabs(n,k)
    for k in range(ngenomes):
        for n in range(ngenomes):
            if n <= ngenomes - k:
                Pabs[n, k] = np.prod((ngenomes - k - np.arange(1, n + 2)) /
                                     (ngenomes + 1 - np.arange(1, n + 2)))

    #calculate Pall(n,k)
    for n in range(ngenomes):
        for k in range(ngenomes):
            Pprod = np.zeros(n + 1)
            for m in range(n + 1):
                Pprod[m] = (k - m + 1) / (ngenomes - m)
            Pall[n, k] = np.prod(Pprod)

    #calculate gk.pan and gk.core
    for n in range(ngenomes):
        gk_pan[n] = np.sum(Gk * (1 - Pabs[n, :]))
        gk_core[n] = np.sum(Gk[n:] * Pall[n, n:])

    return ([gk_pan, gk_core])


def f_pangenome(pa_matrix, nreps):
    pan_mat = np.zeros((pa_matrix.shape[1], nreps))
    for i in range(nreps):
        perm = np.random.choice(pa_matrix.shape[1],
                                size=pa_matrix.shape[1],
                                replace=False)
        pan_mat[:, i] = np.sum(np.cumsum(pa_matrix[:, perm], 1) > 0, 0)

    return (pan_mat)


def f_core(pa_matrix, nreps):
    core_mat = np.zeros((pa_matrix.shape[1], nreps))
    for i in range(nreps):
        perm = np.random.choice(pa_matrix.shape[1],
                                size=pa_matrix.shape[1],
                                replace=False)
        cs = np.cumsum(pa_matrix[:, perm], 1)
        core_mat[:, i] = np.sum(cs == np.arange(1, pa_matrix.shape[1] + 1), 0)

    return (core_mat)


def f_coalescent(ngenomes, rho1, theta1, gess=None, rho2=None, theta2=None):

    if gess is None:
        gess = 0.0
    if rho2 is None:
        rho2 = 1.0
        theta2 = 0.0

    core = np.zeros(ngenomes)
    pangenome = np.zeros(ngenomes)

    for k in range(ngenomes):
        #size of pan-genome after sampling K genomes from N
        pangenome[k] = (theta1 * np.sum(1 / (rho1 - 1 + np.arange(1, k + 2))) +
                        theta2 * np.sum(1 / (rho2 - 1 + np.arange(1, k + 2))) +
                        gess)

        specprod1 = (k - np.arange(1, k + 2) + 2) / (k - np.arange(1, k + 2) +
                                                     rho1 + 1)
        specprod2 = (k - np.arange(1, k + 2) + 2) / (k - np.arange(1, k + 2) +
                                                     rho2 + 1)

        core[k] = ((theta1 / (k + 1)) * np.prod(specprod1[:(k + 1)]) +
                   (theta2 / (k + 1)) * np.prod(specprod2[:(k + 1)]) + gess)

    return ([pangenome, core])


def f_coalescent_spec(ngenomes,
                      rho1,
                      theta1,
                      gess=None,
                      rho2=None,
                      theta2=None):

    if gess is None:
        gess = 0.0
    if rho2 is None:
        rho2 = 1.0
        theta2 = 0.0

    spec = np.zeros(ngenomes)
    specprod1 = (ngenomes - np.arange(1, ngenomes + 1) +
                 1) / (ngenomes - np.arange(1, ngenomes + 1) + rho1)
    if theta2 > 0:
        specprod2 = (ngenomes - np.arange(1, ngenomes + 1) +
                     1) / (ngenomes - np.arange(1, ngenomes + 1) + rho2)
    else:
        specprod2 = np.zeros(ngenomes)

    for k in range(ngenomes):
        spec[k] = ((theta1 / (k + 1)) * np.prod(specprod1[:(k + 1)]) +
                   (theta2 / (k + 1)) * np.prod(specprod2[:(k + 1)]))

    spec[ngenomes - 1] = spec[ngenomes - 1] + gess

    return (spec)


def get_tree_table(tree):
    # generate array form of tree
    # prepare matrix representation of tree
    nnodes = 0
    for node in tree.leaf_node_iter():
        node.label = nnodes
        nnodes += 1
    for node in tree.postorder_internal_node_iter():
        node.label = nnodes
        nnodes += 1
    tree_array = np.zeros((nnodes, 9))
    leaves = []
    node_index = {}
    for i, node in enumerate(tree.leaf_node_iter()):
        leaves.append(i)
        node_index[node.label] = i
        tree_array[i][0] = -1
        tree_array[i][1] = -1
        tree_array[i][3] = node.edge.length
        tree_array[i][4] = 0  # number of descendants

    nleaves = len(leaves)
    for i, node in enumerate(tree.postorder_internal_node_iter()):
        j = i + nleaves
        node_index[node.label] = j
        children = node.child_nodes()

        tree_array[j][0] = node_index[children[0].label]
        tree_array[j][1] = node_index[children[1].label]
        tree_array[j][3] = node.edge.length
        tree_array[j][4] = len(node.leaf_nodes()) + len(
            list(node.postorder_internal_node_iter())) - 1

    return (tree_array)


def f_fixed_spec(treetable, v1, u1, gess=None, v2=None, u2=None):

    if gess is None:
        gess = 0.0
    if v2 is None:
        v2 = 1.0
        u2 = 0.0

    treetable[:, 5] = 1 - np.exp(-v1 * treetable[:, 3])
    treetable[:, 6] = 1 - np.exp(-v2 * treetable[:, 3])
    treetable[:, 7] = np.exp(-v1 * treetable[:, 3])
    treetable[:, 8] = np.exp(-v2 * treetable[:, 3])

    nn = treetable.shape[0]  #number of nodes
    ng = int((nn + 1) / 2)  #number of genomes if bifurcating tree
    nd = np.zeros(nn,
                  dtype=int)  #number of descendent nodes not counting itself
    ndg = np.zeros(nn, dtype=int)  #number of descendent genomes
    gexp1 = np.zeros(nn)  #expected number of genes
    gexp2 = np.zeros(nn)  #expected number of genes
    gprob1 = np.zeros((nn, ng + 1))  #probability of k genomes from k=0
    gprob2 = np.zeros((nn, ng + 1))  #probability of k genomes from k=0

    for a in range(nn):
        nd[a] = treetable[a, 4]
        ndg[a] = (nd[a] + 2) / 2  #number of descendent genomes (tips)

        #if root else
        if a == (nn - 1):
            gexp1[a] = u1 / v1
            gexp2[a] = u2 / v2
        else:
            gexp1[a] = (u1 / v1) * treetable[a, 5]
            gexp2[a] = (u2 / v2) * treetable[a, 6]

        k = 0
        if treetable[a][4] == 0:
            gprob1[a, k + 1] = 1
            gprob2[a, k + 1] = 1
        else:
            desc_a = int(treetable[a, 0])
            desc_b = int(treetable[a, 1])

            #for k = 0
            gprob1[a, k] = ((treetable[desc_a, 5] +
                             treetable[desc_a, 7] * gprob1[desc_a, k]) *
                            (treetable[desc_b, 5] +
                             treetable[desc_b, 7] * gprob1[desc_b, k]))
            gprob2[a, k] = ((treetable[desc_a, 6] +
                             treetable[desc_a, 8] * gprob2[desc_a, k]) *
                            (treetable[desc_b, 6] +
                             treetable[desc_b, 8] * gprob2[desc_b, k]))

            #for 1 <= k <= ndg[a]
            for k in range(ndg[a]):
                pre_sum = 0
                for j in range(0, k + 2):
                    pre_sum += gprob1[desc_a, j] * gprob1[desc_b, k - j + 1]

                gprob1[a, k + 1] = (
                    treetable[desc_a, 7] * treetable[desc_b, 5] *
                    gprob1[desc_a, k + 1] + treetable[desc_b, 7] *
                    treetable[desc_a, 5] * gprob1[desc_b, k + 1] +
                    treetable[desc_a, 7] * treetable[desc_b, 7] * pre_sum)

                pre_sum = 0
                for j in range(0, k + 2):
                    pre_sum += gprob2[desc_a, j] * gprob2[desc_b, k - j + 1]
                gprob2[a, k + 1] = (
                    treetable[desc_a, 8] * treetable[desc_b, 6] *
                    gprob2[desc_a, k + 1] + treetable[desc_b, 8] *
                    treetable[desc_a, 6] * gprob2[desc_b, k + 1] +
                    treetable[desc_a, 8] * treetable[desc_b, 8] * pre_sum)

    Gk1 = np.zeros(ng)
    for k in range(ng):
        for a in range(nn):
            Gk1[k] += gexp1[a] * gprob1[a, k + 1]

    Gk2 = np.zeros(ng)
    for k in range(ng):
        for a in range(nn):
            Gk2[k] += gexp2[a] * gprob2[a, k + 1]

    Gk = Gk1 + Gk2
    Gk[ng - 1] = Gk[ng - 1] + gess

    return (Gk)


def f_theory_dist(params,
                  data,
                  constr,
                  modeltype,
                  fit,
                  genomesize=None,
                  ng=None,
                  treetable=None,
                  return_theory=False):

    if np.sum(params < 0) > 0:
        return (math.inf)

    if constr:
        if len(params) == 1:  # 1D constrained
            v1 = params[0]
            u1 = v1 * genomesize
            v2 = 255
            gess = 0
            u2 = 0
        elif len(params) == 2:  # 1D+E constrained
            v1 = params[0]
            gess = np.power(10, params[1])
            u1 = v1 * (genomesize - gess)
            v2 = 255
            u2 = 0
        elif len(params) == 3:  # 2D constrained
            v1 = params[0]
            u1 = params[1]
            v2 = params[2]
            gess = 0
            u2 = v2 * (genomesize - u1 / v1)
        elif len(params) == 4:  # 2D+E constrained
            v1 = params[0]
            u1 = params[1]
            gess = np.power(10, params[2])
            v2 = params[3]
            u2 = v2 * (genomesize - u1 / v1 - gess)
    else:
        if len(params) == 2:  # 1D unconstrained
            v1 = params[0]
            u1 = params[1]
            v2 = 255
            gess = 0
            u2 = 0
        elif len(params) == 3:  # 1D+E unconstrained
            v1 = params[0]
            u1 = params[1]
            gess = np.power(10, params[2])
            v2 = 255
            u2 = 0
        elif len(params) == 4:  # 2D unconstrained
            v1 = params[0]
            u1 = params[1]
            gess = 0
            v2 = params[2]
            u2 = params[3]
        elif len(params) == 5:  # 2D+E unconstrained
            v1 = params[0]
            u1 = params[1]
            gess = np.power(10, params[2])
            v2 = params[3]
            u2 = params[4]

    if modeltype == "fixed" and fit == "gf":
        theory = f_fixed_spec(treetable, v1, u1, gess, v2, u2)
    elif modeltype == "fixed" and fit == "cp":
        result = f_meanpancore(f_fixed_spec(treetable, v1, u1, gess, v2, u2))
        theory = np.concatenate(result)
    elif modeltype == "coalescent" and fit == "gf":
        theory = f_coalescent_spec(ng, v1, u1, gess, v2, u2)
    elif modeltype == "coalescent" and fit == "cp":
        result = f_coalescent(ng, v1, u1, gess, v2, u2)
        theory = np.concatenate(result)
    else:
        raise RuntimeError("Invalid fit parameters!")

    if return_theory:
        return (result)

    return (np.sqrt(np.sum(np.power(theory - data, 2)) / len(theory)))


def get_options():
    import argparse

    description = 'Estimate model parameters for either the Infinitely Many Genes Model using gene frequencies'
    parser = argparse.ArgumentParser(description=description,
                                     prog='panaroo_img_est')

    parser.add_argument("--tree",
                        dest="tree",
                        help="A dated phylogeny.",
                        type=str)

    parser.add_argument("--pa",
                        dest="presence_absence",
                        help="A presence/absence produced by Panaroo.",
                        type=str,
                        required=True)

    parser.add_argument("-o",
                        "--out_dir",
                        dest="output_dir",
                        required=True,
                        help="location of an output directory",
                        type=lambda x: is_valid_folder(parser, x))

    parser.add_argument("-D",
                        dest="n_classes",
                        help=("Number of seperate rate classes to use for " +
                              "the dispensable genome. Can be either 1 or 2."),
                        type=int,
                        choices={1, 2},
                        required=True)

    parser.add_argument("--no_essential",
                        dest="no_essential",
                        help=("Removes essential gene class from model"),
                        action='store_true',
                        default=False)

    parser.add_argument(
        "--no_constraint",
        dest="no_constraint",
        help=("Removes constraint that u/v must equal the genome size."),
        action='store_true',
        default=False)

    parser.add_argument(
        "--model",
        dest="model",
        help=("Model to fit. Can be one of 'coalescent' or 'fixed'."),
        type=str,
        choices={'coalescent', 'fixed'},
        default='fixed')

    parser.add_argument(
        "--fit",
        dest="fit",
        help=("Whether to use the gene frequency spectrum or the" +
              " core/pangeome size for fitting (default=gf)"),
        type=str,
        choices={'gf', 'cp'},
        default='gf')

    parser.add_argument("--init_u1",
                        dest="u1",
                        help=("Initial value for u1 (default=0.01)."),
                        type=float,
                        default=0.01)

    parser.add_argument("--init_u2",
                        dest="u2",
                        help=("Initial value for u2 (default=0.01)."),
                        type=float,
                        default=0.01)

    parser.add_argument("--init_v1",
                        dest="v1",
                        help=("Initial value for v1 (default=0.01)."),
                        type=float,
                        default=0.01)

    parser.add_argument("--init_v2",
                        dest="v2",
                        help=("Initial value for v2 (default=0.01)."),
                        type=float,
                        default=0.01)

    parser.add_argument(
        "--init_ess",
        dest="gess",
        help=(
            "Initial value for the number of essential genes (default=2000)."),
        type=float,
        default=2000)

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

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # load input data
    isolates, genes, pa_matrix = load_pa(args.presence_absence)
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

    genomesize = np.mean(np.sum(pa_matrix > 0, 0))  # mean genome size
    ngenomes = pa_matrix.shape[1]

    Gk = f_getspectrum(pa_matrix)
    if args.fit == "cp":
        data = np.concatenate(f_meanpancore(Gk))
    else:
        data = Gk

    if args.no_constraint:
        params_init = [args.v1, args.u1]
        bounds = [(1e-7, 1e3), (1e-7, 1e3)]
        if not args.no_essential:
            params_init += [np.log10(args.gess)]
            bounds += [(1e-7, 7)]
        if args.n_classes == 2:
            params_init += [args.v2, args.u2]
            bounds += [(1e-7, 1e3), (1e-7, 1e3)]
    else:
        params_init = [args.v1]
        bounds = [(1e-7, 1e3)]
        if args.n_classes == 2:
            params_init += [args.u1]
            bounds += [(1e-7, 1e3)]
        if not args.no_essential:
            params_init += [np.log10(args.gess)]
            bounds += [(1e-7, 7)]
        if args.n_classes == 2:
            params_init += [args.v2]
            bounds += [(1e-7, 1e3)]

    options = {"maxiter": 1e4}

    if args.model == 'coalescent':
        # result = optimize.minimize(f_theory_dist,
        #                         bounds=bounds,
        #                         x0=params_init,
        #                         method='L-BFGS-B',
        #                         args=(data, (not args.no_constraint),
        #                                 args.model, args.fit, genomesize,
        #                                 ngenomes),
        #                         options=options)

        result = optimize.basinhopping(f_theory_dist,
                                       x0=params_init,
                                       niter=10,
                                       minimizer_kwargs={
                                           'args':
                                           (data, (not args.no_constraint),
                                            args.model, args.fit, genomesize,
                                            ngenomes),
                                           "method":
                                           'L-BFGS-B',
                                           "bounds":
                                           bounds
                                       })

    else:
        if args.tree is None:
            raise RuntimeError("A phylogeny is required for the fixed model!")
        tree_table = get_tree_table(tree)

        # result = optimize.minimize(f_theory_dist,
        #                         bounds=bounds,
        #                         x0=params_init,
        #                         method='L-BFGS-B',
        #                         args=(data, (not args.no_constraint),
        #                                 args.model, args.fit, genomesize,
        #                                 ngenomes, tree_table),
        #                         options=options)

        result = optimize.basinhopping(f_theory_dist,
                                       x0=params_init,
                                       niter=10,
                                       minimizer_kwargs={
                                           'args':
                                           (data, (not args.no_constraint),
                                            args.model, args.fit, genomesize,
                                            ngenomes, tree_table),
                                           "method":
                                           'L-BFGS-B',
                                           "bounds":
                                           bounds
                                       })

    if args.verbose:
        print("***** Optimisation Summary ******")
        print(result)
        print("*********************************\n\n")

    # summarise results
    params = list(result['x'])
    if args.no_constraint:
        result_summary = [('v1', params.pop(0)), ('u1', params.pop(0))]
        if not args.no_essential:
            result_summary += [('ess', np.power(10, params.pop(0)))]
        if args.n_classes == 2:
            result_summary += [('v2', params.pop(0)), ('u2', params.pop(0))]
    else:
        result_summary = [('v1', params.pop(0))]
        if args.n_classes == 2:
            result_summary += [('u1', params.pop(0))]
        if not args.no_essential:
            result_summary += [('ess', np.power(10, params.pop(0)))]
        if args.n_classes == 2:
            result_summary += [('v2', params.pop(0))]

    if args.verbose:
        print("***** Inferred Parameters ******")
        for p in result_summary:
            print(p[0] + ": " + str(p[1]))
        print("*********************************\n\n")

    with open(args.output_dir + "infered_parameters.csv", 'w') as outfile:
        outfile.write("parameter,value\n")
        for p in result_summary:
            outfile.write(p[0] + "," + str(p[1]) + "\n")

    pan, core = f_theory_dist(result['x'],
                              data, (not args.no_constraint),
                              args.model,
                              "cp",
                              genomesize=genomesize,
                              ng=ngenomes,
                              treetable=tree_table,
                              return_theory=True)

    with open(args.output_dir + "core_pangenome_size_estimate.csv",
              'w') as outfile:
        outfile.write("genome,pan,core\n")
        for i, p, c in zip(range(len(pan)), pan, core):
            outfile.write(str(i) + "," + str(p) + "," + str(c) + "\n")

    return


if __name__ == '__main__':
    main()
