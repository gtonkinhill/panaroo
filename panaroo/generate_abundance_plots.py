import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

import os
import argparse
from .isvalid import *
from .__init__ import __version__


def read_presence_absence(filename):
    # load matrix into binary array
    matrix_txt = np.loadtxt(filename, dtype=str, delimiter=",", comments=None)
    sample_names = matrix_txt[0, 15:]
    gene_names = matrix_txt[0, 1:]
    pa_matrix = matrix_txt[1:, 15:] != ""
    return (pa_matrix, gene_names, sample_names)


def get_curve_w_ci(pa_matrix, n_boot=100, method="chao2"):

    n_samples = np.size(pa_matrix, 1)

    stat_ci = []

    if method == "ICE":
        start = 20
    else:
        start = 2

    for n in range(start, n_samples + 1):

        if method == 'acc':
            qm_boot = [
                pa_matrix[:,
                          np.random.choice(n_samples, size=n, replace=True)]
                for i in range(1, n_boot + 1)
            ]
        else:
            qm_boot = [
                get_q_m(pa_matrix[:,
                                  np.random.
                                  choice(n_samples, size=n, replace=True)])
                for i in range(1, n_boot + 1)
            ]

        if method == "chao2":
            stat_boot = [chao2(qm[0], qm[1], qm[2]) for qm in qm_boot]
        elif method == "ICE":
            stat_boot = [ICE(qm[0], pa_matrix, n_samples) for qm in qm_boot]
        elif method == "jack1":
            stat_boot = [jackknife(qm[0], qm[1], qm[2], 1) for qm in qm_boot]
        elif method == "jack2":
            stat_boot = [jackknife(qm[0], qm[1], qm[2], 2) for qm in qm_boot]
        else:
            stat_boot = [acc_curve(qm) for qm in qm_boot]

        Q = np.quantile(stat_boot, [0.025, 0.5, 0.975])
        Q[1] = np.mean(stat_boot)
        stat_ci.append(Q)

    return range(start, n_samples + 1), stat_ci


def get_q_m(pa_matrix):
    # q_k be the number of genes present in exactly k genomes
    f = np.sum(pa_matrix, axis=1)
    q = np.bincount(f)
    m = 0.0
    sobs = np.sum(f > 0)
    for i in range(1, np.size(q)):
        m += q[i]
    return q, m, sobs


def chao2(q, m, sobs):
    # calculates the Chao 2 statistic (for replicated incidence data)
    c2 = sobs + ((m - 1) / m) * (q[1] * (q[1] - 1)) / (2 * (q[2] + 1))
    return c2


def acc_curve(pa_matrix):
    # calculates the number of unique genes
    return np.sum(np.sum(pa_matrix, axis=1) > 0)


def ICE(q, pa_matrix, n_samples, thresh=10):
    # ICE estimator of species richness
    S_infr = np.sum(q[1:(thresh + 1)])
    S_freq = np.sum(q[(thresh + 1):(n_samples + 1)])
    C_infr = 1 - q[1] / sum([i * q[i] for i in range(1, (thresh + 1))])
    f = np.sum(pa_matrix, axis=1)
    T_infr = np.sum(np.sum(pa_matrix[f <= 10, :], axis=0) > 0)
    iQ = sum([i * q[i] for i in range(1, (thresh + 1))])
    i_iQ = sum([i * (i - 1) * q[i] - 1 for i in range(1, (thresh + 1))])
    lambda_ice = max((S_infr / C_infr) * (T_infr / (T_infr - 1)) *
                     (i_iQ / (iQ * (iQ - 1))) - 1, 0)

    return S_freq + S_infr / C_infr + q[1] / C_infr * lambda_ice


def jackknife(q, m, sobs, order=1):
    # first and second-order jackknife richness estimator
    if order == 1:
        S = sobs + q[1] * ((m - 1) / m)
    else:
        S = sobs + ((q[1] * (2 * m - 3)) / m - ((q[2] * np.power(m - 2, 2)) /
                                                (m * (m - 1))))

    return S


def plot_quantiles(quantiles,
                   sample_nums,
                   method,
                   outdir,
                   color_shading='#1f78b4',
                   color_mean='#1f78b4'):

    plt.style.use('ggplot')

    fig = plt.figure()

    plt.fill_between(sample_nums, [q[0] for q in quantiles],
                     [q[2] for q in quantiles],
                     color=color_shading,
                     alpha=0.2)
    # plot the mean on top
    plt.plot(sample_nums, [q[1] for q in quantiles], color_mean)

    plt.grid(True)
    plt.xlabel("Number of Samples")
    plt.ylabel(method + " Diversity")
    plt.tight_layout()

    fig.savefig(outdir + method + ".png")

    with open(outdir + method + "_data.csv", 'w') as outfile:
        outfile.write("0.025,median,0.975\n")
        for quant in quantiles:
            outfile.write(",".join([str(q) for q in quant]) + "\n")

    return


def generate_plot(pa_file, method, nboot, outdir):

    pa_matrix, gene_names, sample_names = read_presence_absence(pa_file)

    sample_nums, quantiles = get_curve_w_ci(pa_matrix,
                                            n_boot=nboot,
                                            method=method)

    plot_quantiles(quantiles, sample_nums, method, outdir)

    return


def get_options():
    import argparse

    description = 'Generate summary plots from Panaroo run'
    parser = argparse.ArgumentParser(description=description,
                                     prog='plot_panaroo_summary')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-i",
        "--input",
        dest="pa_file",
        required=True,
        help="gene presence absence file generated by Panaroo")
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=lambda x: is_valid_folder(parser, x))

    # Other options
    parser.add_argument(
        "--graph_type",
        dest="graph_type",
        help="the type of graph to generate",
        choices={'all', 'chao2', 'ICE', 'jack1', 'jack2', 'acc'})
    parser.add_argument(
        "--nboot",
        dest="nboot",
        type=int,
        default=100,
        help="the number of bootstrap replicates (default=100)")
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    if args.graph_type == 'all':
        for method in ['chao2', 'ICE', 'jack1', 'jack2', 'acc']:
            generate_plot(args.pa_file, method, args.nboot, args.output_dir)
    else:
        generate_plot(args.pa_file, args.graph_type, args.nboot,
                      args.output_dir)

    return


if __name__ == '__main__':
    main()
