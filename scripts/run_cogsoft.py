import argparse
import subprocess
import sys, os
import shutil
import gffutils as gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import re, string
import numpy as np
from collections import defaultdict

from panaroo.prokka import get_gene_sequences

def prepare_gffs(gff_file_names, outdir):
    
    out_file = outdir+"all_proteins.fasta"

    genomeids = []

    mapping = {}
    seen = set()

    with open(outdir + "all_proteins.fasta", 'w') as outfile:
        for i, gfffile in enumerate(gff_file_names):
            prefix = os.path.splitext(os.path.basename(gfffile))[0]
            genomeids.append(re.sub('[\W_]+', '', prefix))
            gff_file = open(gfffile, 'r')
            seq_dict, prot_seq_dict = get_gene_sequences(gff_file, i)
            SeqIO.write(prot_seq_dict, outfile, 'fasta')

            for protien in prot_seq_dict:
                clustering_id = protien.id
                relevant_seqrecord = seq_dict[clustering_id]
                if relevant_seqrecord.id in seen:
                    raise RuntimeError("Duplicate gene names!")
                else:
                    seen.add(relevant_seqrecord.id)
                mapping[clustering_id] = relevant_seqrecord.id

    with open(outdir + "genome_to_protein.csv", 'w') as outfile:
        for clustering_id in mapping:
            genome = genomeids[int(clustering_id.split("_")[0])]
            outfile.write(clustering_id + "," + genome
                 + "\n")

    return mapping, genomeids


def run_blast(outdir, ncpus, verbose=True):

    # format a BLASTable database
    cmd = "makeblastdb "
    cmd += " -in " + outdir + "all_proteins.fasta" 
    cmd += " -dbtype prot"
    cmd += " -out " + outdir + "BLASTable" 

    if verbose:
        print("running cmd: ", cmd)
    subprocess.run(cmd, shell=True, check=True)

    # unfiltered BLAST results in the ./BLASTno/ directory
    cmd = "psiblast"
    cmd += " -num_threads " + str(ncpus)
    cmd += " -query " + outdir + "all_proteins.fasta" 
    cmd += " -db " + outdir + "BLASTable" 
    cmd += " -show_gis -outfmt 7 -num_descriptions 1000"
    cmd += " -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no"
    cmd += " -out " + outdir + "/BLASTno/ThreeByThree.tab"

    if verbose:
        print("running cmd: ", cmd)
    if not os.path.exists(outdir + "BLASTno"):
        os.mkdir(outdir + "BLASTno")
    subprocess.run(cmd, shell=True, check=True)

    # filtered BLAST results in the ./BLASTff/ directory
    cmd = "psiblast"
    cmd += " -num_threads " + str(ncpus)
    cmd += " -query " + outdir + "all_proteins.fasta"
    cmd += " -db " + outdir + "BLASTable" 
    cmd += " -show_gis -outfmt 7 -num_descriptions 1000"
    cmd += " -num_alignments 1000 -dbsize 100000000" 
    cmd += " -comp_based_stats T -seg yes"
    cmd += " -out " + outdir + "BLASTff/ThreeByThree.tab "

    if verbose:
        print("running cmd: ", cmd)
    if not os.path.exists(outdir + "BLASTff"):
        os.mkdir(outdir + "BLASTff")
    subprocess.run(cmd, shell=True, check=True)

    return

def run_diamond(outdir, ncpus, verbose=True):

    # format a DIAMOND data base
    cmd = "diamond makedb"
    cmd += " --in " + outdir + "all_proteins.fasta"
    cmd += " --db " + outdir + "DIAMONDTable" 

    if verbose:
        print("running cmd: ", cmd)
    subprocess.run(cmd, shell=True, check=True)

    # unfiltered DIAMOND
    cmd = "diamond blastp"
    cmd += " --threads " + str(ncpus)
    cmd += " --query " + outdir + "all_proteins.fasta" 
    cmd += " --db " + outdir + "DIAMONDTable" 
    cmd += " --outfmt 6"
    cmd += " --max-target-seqs 1000"
    cmd += " --masking 0"
    cmd += " --out " + outdir + "/BLASTno/ThreeByThree.tab"

    if verbose:
        print("running cmd: ", cmd)
    if not os.path.exists(outdir + "BLASTno"):
        os.mkdir(outdir + "BLASTno")
    subprocess.run(cmd, shell=True, check=True)

    # filtered DIAMOND results in the ./BLASTff/ directory
    cmd = "diamond blastp"
    cmd += " --threads " + str(ncpus)
    cmd += " --query " + outdir + "all_proteins.fasta" 
    cmd += " --db " + outdir + "DIAMONDTable" 
    cmd += " --outfmt 6"
    cmd += " --max-target-seqs 1000"
    cmd += " --masking 1"
    cmd += " --out " + outdir + "BLASTff/ThreeByThree.tab "

    if verbose:
        print("running cmd: ", cmd)
    if not os.path.exists(outdir + "BLASTff"):
        os.mkdir(outdir + "BLASTff")
    subprocess.run(cmd, shell=True, check=True)

    return

def run_cog_soft(outdir, genomeids, verbose=True):

    # Preparation of the "sequence Universe"
    # i.e make the ./BLASTconv/hash.csv file 
    cmd = "COGmakehash"
    cmd += " -i=" + outdir + "genome_to_protein.csv"
    cmd += " -o=" + outdir + "BLASTconv"
    cmd += " -s=',' -n=1"

    if not os.path.exists(outdir + "BLASTconv"):
        os.mkdir(outdir + "BLASTconv")
    if verbose:
        print("running cmd: ", cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Processing of BLAST results
    cmd = "COGreadblast"
    cmd += " -d=" + outdir + "BLASTconv"
    cmd += " -u=" + outdir + "BLASTno"
    cmd += " -f=" + outdir + "BLASTff"
    cmd += " -s=" + outdir + "BLASTno"
    cmd += " -e=0.1 -q=2 -t=2"

    if verbose:
        print("running cmd: ", cmd)
    subprocess.run(cmd, shell=True, check=True)


    # Lineage specific expansions

    # here we use the BLAST results to choose a rough outgroup for each genome
    nsamples = len(genomeids)
    hit_counts = np.zeros((nsamples, nsamples), dtype=int)
    with open(outdir + "BLASTff/ThreeByThree.tab", 'r') as infile:
        for line in infile:
            if line[0]=="#": continue
            line = line.strip().split()
            hitA = int(line[0].split("_")[0])
            hitB = int(line[1].split("_")[0])
            hit_counts[hitA, hitB] += 1
            hit_counts[hitB, hitA] += 1

    # now create the required outgroup file
    with open(outdir + "GenThree.job.csv", 'w') as outfile:
        for i in range(nsamples):
            outfile.write(genomeids[i] + "," + genomeids[np.argmin(hit_counts[i,:])] + "\n")

    # now identify lineage specific expansions
    cur_dir = os.getcwd()
    os.chdir(outdir)

    cmd = "COGlse" 
    cmd += " -d=" + "BLASTconv"
    cmd += " -j=" + "GenThree.job.csv"
    cmd += " -p=" + "genome_to_protein.csv"
    cmd += " -o=" + "GenThree.lse.csv"

    if not os.path.exists("data"):
        os.mkdir("data")
    if not os.path.exists("data/tmp"):
        os.mkdir("data/tmp")
    if verbose:
        print("running cmd: ", cmd)
    subprocess.run(cmd, shell=True, check=True)

    os.chdir(cur_dir)

    # Making clusters from symmetrical best hits
    cmd = "COGtriangles"
    cmd += " -i=" + outdir + "BLASTconv"
    cmd += " -q=" + outdir + "genome_to_protein.csv"
    cmd += " -l=" + outdir + "GenThree.lse.csv"
    cmd += " -o=" + outdir + "GenThree.cls.csv"
    cmd += " -t=0.5 -e=0.01 -n='CLS' -s=1"

    if verbose:
        print("running cmd: ", cmd)
    subprocess.run(cmd, shell=True, check=True)

    return



def post_process_fmt(outdir, gff_file_names, mapping):

    # get original filenames
    prefixes = [os.path.splitext(os.path.basename(gfffile))[0] for gfffile in gff_file_names]

    # generate an index for files
    n_samples = len(prefixes)
    
    # collect clusters
    clusters = defaultdict(list)
    with open(outdir + "GenThree.cls.csv", 'r') as infile:
        for line in infile:
            line=line.strip().split(",")
            clusters[line[6]].append(line[0])

    # write out in roary like format
    with open(outdir + "COGsoft_gene_presence_absence.csv", 'w') as outfile:
        outfile.write(",".join(prefixes) + "\n")
        for cluster_name in clusters:
            pa = n_samples*[""]
            cluster = clusters[cluster_name]
            for g in cluster:
                pa[int(g.split("_")[0])] = mapping[g]
            outfile.write(",".join(pa) + "\n")

    return



def main():
    parser = argparse.ArgumentParser(description="""Runs COGsoft on GFF3 files.""")
    parser.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=str)

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="input GFF3 files (usually output from running Prokka)",
        type=str,
        nargs='+')

    parser.add_argument(
        "--aligner",
        dest="alr",
        help=
        "Specify an aligner. Options:'BLAST', 'DIAMOND' (default)",
        type=str,
        choices=['DIAMOND', 'BLAST'],
        default="DIAMOND")

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    args = parser.parse_args()

    args.output_dir = os.path.join(args.output_dir, "")

    mapping, genomeids = prepare_gffs(args.input_files, args.output_dir)

    if args.alr=="BLAST":
        run_blast(args.output_dir, args.n_cpu)
    else:
        run_diamond(args.output_dir, args.n_cpu)

    run_cog_soft(args.output_dir, genomeids)

    post_process_fmt(args.output_dir, args.input_files, mapping)

    # clean up BLAST/DIAMOND output as it uses too much space
    shutil.rmtree(args.output_dir + "BLASTff")
    shutil.rmtree(args.output_dir + "BLASTno")
    shutil.rmtree(args.output_dir + "BLASTconv")

    return

if __name__ == '__main__':
    main()