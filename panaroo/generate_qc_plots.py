import matplotlib
from matplotlib import pyplot as plt
import subprocess
import tempfile
import os, sys
from joblib import Parallel, delayed
import numpy as np
from sklearn import manifold
from collections import Counter
import plotly.graph_objs as go
import plotly.offline as offline
from plotly import tools

from .isvalid import *
from .__init__ import __version__


def get_mash_dist(input_gffs, outdir, n_cpu=1, quiet=True):

    # build mash sketch
    mash_cmd = "mash triangle"
    mash_cmd += " -p " + str(n_cpu)

    #Set up two lists of gffs to input into mash. This is a bit messy but
    # allows for an arbitary number of gffs.
    temp_output_file1 = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_output_file2 = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_output_file1.close()
    temp_output_file2.close()
    
    with open(temp_output_file1.name, 'w') as outfile:
        outfile.write(input_gffs[0])
    with open(temp_output_file2.name, 'w') as outfile:
        for gff in input_gffs[1:]:
            outfile.write(gff + "\n")
    
    mash_cmd += " -l " + temp_output_file1.name
    mash_cmd += " " + temp_output_file2.name
    mash_cmd += " > " + outdir + "mash_dist.txt"

    if not quiet:
        print("running cmd: " + mash_cmd)

    subprocess.run(mash_cmd, shell=True, check=True)

    # load distance matrix
    dist_mat = np.zeros((len(input_gffs), len(input_gffs)), dtype=float)
    with open(outdir + "mash_dist.txt", 'r') as infile:
        next(infile)
        for i, line in enumerate(infile):
            line = line.strip().split()
            for j, d in enumerate(line[1:]):
                dist_mat[i, j] = float(d)
                dist_mat[j, i] = dist_mat[i, j]

    # get simplified file names
    file_names = [
        os.path.splitext(os.path.basename(gff))[0] for gff in input_gffs
    ]

    # clean up
    os.remove(temp_output_file1.name)
    os.remove(temp_output_file2.name)

    return dist_mat, file_names


def plot_MDS(dist_mat, file_names, outdir, no_plot=False):

    # get MDS projection
    mds = manifold.MDS(n_components=2, dissimilarity="precomputed")
    projection = mds.fit(dist_mat)
    coords = projection.embedding_

    #write MDS coordinates to disk
    with open(outdir + "mds_coords.txt", "w") as contig_out:
        contig_out.write("sample\tcoordx\tcoordy\n")
        for i, coord in zip(file_names, coords):
            contig_out.write("%s\t%s\t%s\n" % (i, coord[0], coord[1]))

    if no_plot: return

    # find margins for plot
    c_min = np.min(coords) - abs(np.quantile(coords, 0.05))
    c_max = np.max(coords) + abs(np.quantile(coords, 0.05))

    # generate static plot
    plt.style.use('ggplot')
    fig = plt.figure()
    plt.scatter(coords[:, 0], coords[:, 1])
    plt.grid(True)
    plt.xlabel("MDS Dimension 1")
    plt.ylabel("MDS Dimension 2")
    plt.xlim((c_min, c_max))
    plt.ylim((c_min, c_max))
    plt.tight_layout()
    fig.savefig(outdir + "MDS_mash_plot.png")

    # generate interactive plot
    trace = go.Scatter(x=coords[:, 0],
                       y=coords[:, 1],
                       text=file_names,
                       mode='markers')
    layout = go.Layout(xaxis=dict(autorange=True,
                                  showgrid=True,
                                  zeroline=True,
                                  showline=False,
                                  ticks='',
                                  range=[c_min, c_max],
                                  type="linear",
                                  exponentformat="SI",
                                  showexponent='none',
                                  showticklabels=True),
                       yaxis=dict(autorange=True,
                                  showgrid=True,
                                  zeroline=True,
                                  showline=False,
                                  ticks='',
                                  range=[c_min, c_max],
                                  type="linear",
                                  exponentformat="SI",
                                  showexponent='none',
                                  showticklabels=True))
    data = [trace]
    fig = go.Figure(data=data, layout=layout)
    offline.plot(fig, filename=outdir + "MDS_mash_plot.html", auto_open=False)

    return


def plot_ngenes(input_gffs, outdir, no_plot=True):

    # get simplified file names
    file_names = [
        os.path.splitext(os.path.basename(gff))[0] for gff in input_gffs
    ]

    # count genes
    ngenes = np.zeros(len(input_gffs))
    for i, gff_file in enumerate(input_gffs):
        with open(gff_file, 'r') as gff:
            for line in gff:
                if "##FASTA" in line: break
                if "##" == line[:2]: continue
                if "CDS" not in line: continue
                ngenes[i] += 1

    with open(outdir + "ngenes.txt", "w") as genes_out:
        genes_out.write("sample\tno_genes\n")
        for i, j in zip(file_names, ngenes):
            genes_out.write("%s\t%s\n" % (i, j))

    if no_plot: return

    # generate static plot
    plt.style.use('ggplot')
    fig = plt.figure()
    plt.barh(np.arange(len(ngenes)), ngenes)
    plt.yticks(np.arange(len(ngenes)), file_names)
    plt.grid(True)
    plt.ylabel("File Name")
    plt.xlabel("Number of Genes")
    plt.tight_layout()
    fig.savefig(outdir + "ngenes_barplot.png")

    # generate interactive boxplot
    data = [
        go.Box(y=ngenes,
               text=file_names,
               hoverinfo="text",
               boxpoints='all',
               jitter=0.3,
               pointpos=-1.8)
    ]
    layout = go.Layout(autosize=True,
                       xaxis=dict(title=dict(text='', font=dict(size=18, color='black')),
                                  showticklabels=False,
                                  automargin=True),
                       yaxis=dict(title=dict(text="Number of Genes", font=dict(size=18, color='black')),
                                  showticklabels=True,
                                  tickfont=dict(size=10, color='black')))

    fig = go.Figure(data=data, layout=layout)
    offline.plot(fig, filename=outdir + "ngenes_boxplot.html", auto_open=False)

    return


def plot_ncontigs(input_gffs, outdir, no_plot=False):

    # get simplified file names
    file_names = [
        os.path.splitext(os.path.basename(gff))[0] for gff in input_gffs
    ]

    # count genes
    ncontigs = np.zeros(len(input_gffs))
    for i, gff_file in enumerate(input_gffs):
        with open(gff_file, 'r') as gff:
            in_fasta = False
            for line in gff:
                if in_fasta and (line[0] == ">"):
                    ncontigs[i] += 1
                if "##FASTA" in line:
                    in_fasta = True

    # generate static plot
    with open(outdir + "ncontigs.txt", "w") as contig_out:
        contig_out.write("sample\tno_contigs\n")
        for i, j in zip(file_names, ncontigs):
            contig_out.write("%s\t%s\n" % (i, j))

    if no_plot: return
    
    plt.style.use('ggplot')
    fig = plt.figure()
    plt.barh(np.arange(len(ncontigs)), ncontigs)
    plt.yticks(np.arange(len(ncontigs)), file_names)
    plt.grid(True)
    plt.ylabel("File Name")
    plt.xlabel("Number of Contigs")
    plt.tight_layout()
    fig.savefig(outdir + "ncontigs_barplot.png")

    # generate interactive boxplot
    data = [
        go.Box(y=ncontigs,
               text=file_names,
               hoverinfo="text",
               boxpoints='all',
               jitter=0.3,
               pointpos=-1.8)
    ]
    layout = go.Layout(autosize=True,
                       xaxis=dict(title=dict(text='', 
                                             font=dict(size=18, color='black')),
                                  showticklabels=False,
                                  automargin=True),
                       yaxis=dict(title=dict(text="Number of Contigs", 
                                             font=dict(size=18, color='black')),
                                  showticklabels=True,
                                  tickfont=dict(size=10, color='black')))

    fig = go.Figure(data=data, layout=layout)
    offline.plot(fig,
                 filename=outdir + "ncontigs_boxplot.html",
                 auto_open=False)

    return


def run_mash_screen(gff, mash_ref, outdir):
    # run mash
    mash_cmd = "mash screen -i 0.8 -w"
    mash_cmd += " " + mash_ref

    temp_output_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_output_file.close()
    temp_mash_cmd = mash_cmd + " " + gff + " > " + temp_output_file.name
    subprocess.run(temp_mash_cmd, shell=True, check=True)

    # load output
    hits = []
    with open(temp_output_file.name, 'r') as infile:
        for line in infile:
            line = line.strip().split("\t")
            line[1] = int(line[1].split("/")[0])
            hits.append(tuple(line))
        hits = sorted(hits, key=lambda x: -x[1])[:10]

    # clean up
    os.remove(temp_output_file.name)

    return (hits)


def get_mash_contam(input_gffs, mash_ref, n_cpu, outdir):
    file_names = [
        os.path.splitext(os.path.basename(gff))[0] for gff in input_gffs
    ]

    genome_hits = Parallel(n_jobs=n_cpu)(
        delayed(run_mash_screen)(gff, mash_ref, outdir) for gff in input_gffs)

    with open(outdir + "mash_contamination_hits.tab", 'w') as outfile:
        for i, genome in enumerate(genome_hits):
            for hit in genome:
                outfile.write("\t".join([file_names[i]] +
                                        [str(h) for h in hit]) + "\n")

    return (outdir + "mash_contamination_hits.tab")


def plot_mash_contam(mash_contam_file, outdir):
    # load mash hits from file
    y = []
    x_label = []
    text = []
    hit_count = Counter()
    with open(mash_contam_file, 'r') as infile:
        for line in infile:
            line = line.strip().split("\t")
            y.append(int(line[2]) / 1000.0)
            x_label.append(line[5])
            text.append("File: " + line[0] + "<br>Hit: " + line[6])
            hit_count[line[5]] += int(line[2])

    # sort hits based on prevalence
    hit_index = {}
    tick_labels = []
    tickvals = []
    for i, hit in enumerate(hit_count.most_common()):
        hit_index[hit[0]] = i + 1
        tick_labels.append(hit[0])
        tickvals.append(i + 1)
    x = []
    for xl in x_label:
        jitter = np.random.normal(loc=0, scale=0.1)
        while abs(jitter) > 0.3:
            jitter = np.random.normal(loc=0, scale=0.1)
        x.append(hit_index[xl] + jitter)

    trace = go.Scatter(x=x, y=y, mode='markers', text=text, hoverinfo="text")

    layout = go.Layout(autosize=True,
                       xaxis=dict(title=dict(text='Match', 
                                             font=dict(size=18, color='black')),
                                  showticklabels=True,
                                  tickangle=45,
                                  ticktext=tick_labels,
                                  tickvals=tickvals,
                                  automargin=True,
                                  tickfont=dict(size=8, color='black')),
                       yaxis=dict(title=dict(text="Percentage of shared hash's", 
                                             font=dict(size=18, color='black')),
                                  showticklabels=True,
                                  tickangle=45,
                                  tickfont=dict(size=10, color='black')))

    data = [trace]
    fig = go.Figure(data=data, layout=layout)
    offline.plot(fig,
                 filename=outdir + "mash_contamination_barplot.html",
                 auto_open=False)

    return


def generate_qc_plot(method, input_files, outdir, n_cpu, ref_db=None, no_plot=False):

    # plot MDS
    if method in ["mds", "all"]:
        dist_mat, file_names = get_mash_dist(input_gffs=input_files,
                                             outdir=outdir,
                                             n_cpu=n_cpu,
                                             quiet=True)
        plot_MDS(dist_mat, file_names, outdir, no_plot)

    # plot number of genes
    if method in ["ngenes", "all"]:
        plot_ngenes(input_gffs=input_files, outdir=outdir, no_plot=no_plot)

    # plot number of contigs
    if method in ["ncontigs", "all"]:
        plot_ncontigs(input_gffs=input_files, outdir=outdir, no_plot=no_plot)

    # plot contamination scatter plot
    if (method in ["contam", "all"]):
        if ref_db is None:
            print(
                "No reference mash database given! Skipping contamination plot..."
            )
            print(("One can be downloaded from https://mash.readthedocs.io" +
                   "/en/latest/tutorials.html#screening-a-read-set-for" +
                   "-containment-of-refseq-genomes"))
        else:
            mash_contam_file = get_mash_contam(input_gffs=input_files,
                                               mash_ref=ref_db,
                                               n_cpu=n_cpu,
                                               outdir=outdir)
            if not no_plot:
                plot_mash_contam(mash_contam_file=mash_contam_file, outdir=outdir)

    return


def get_options(args):
    import argparse

    description = 'Generate quality control plots prior to a Panaroo run'
    parser = argparse.ArgumentParser(description=description,
                                     prog='plot_panaroo_qc')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="input GFF3 files (usually output from running Prokka)",
        type=str,
        nargs='+')
    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=lambda x: is_valid_folder(parser, x))

    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--graph_type",
                        dest="graph_type",
                        help="the type of graph to generate (default='all')",
                        choices={'all', 'mds', 'ngenes', 'ncontigs', 'contam'},
                        default="all")
    parser.add_argument("--no_plot",
                        dest="no_plot",
                        help="don't generate the plots. Will only create the data tables",
                        action='store_true',
                        default=False)
    parser.add_argument(
        "--ref_db",
        dest="ref_db",
        help="reference mash database for contamination quantification.")
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args(args)
    return (args)


def main():
    args = get_options(sys.argv[1:])

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    if len(args.input_files) == 1:
        with open(args.input_files[0]) as inhandle:
            args.input_files = inhandle.read().splitlines()
    
    generate_qc_plot(method=args.graph_type,
                     input_files=args.input_files,
                     outdir=args.output_dir,
                     n_cpu=args.n_cpu,
                     ref_db=args.ref_db,
                     no_plot=args.no_plot)

    return


if __name__ == '__main__':
    main()
