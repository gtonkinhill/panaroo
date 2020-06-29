import os
import shutil
import tempfile
import argparse

import mappy as mp
import pyfastx
from collections import defaultdict, Counter
import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from itertools import zip_longest, chain


# def run_minimap2(r1, r2, db, outdir, ncpu=1, quiet=False):
#     # run minimap2 and output to sorted bam
#     cmd = 'minimap2 -ax sr --secondary=yes -N 100000 -m 72'
#     cmd += ' -t ' + str(ncpu) 
#     cmd += ' ' + db
#     cmd += ' ' + r1
#     cmd += ' ' + r2
#     cmd += ' | samtools sort'
#     cmd += ' -@ ' + str(ncpu) 
#     cmd += ' > ' + outdir + 'alignment.bam'

#     if not quiet: print("running cmd: ", cmd)
#     subprocess.run(cmd, shell=True, check=True)

#     # index bam file
#     cmd = 'samtools index ' + outdir + 'alignment.bam'
#     if not quiet: print("running cmd: ", cmd)
#     subprocess.run(cmd, shell=True, check=True)

#     # index db
#     cmd = 'samtools faidx ' + db
#     if not quiet: print("running cmd: ", cmd)
#     subprocess.run(cmd, shell=True, check=True)
#     return(outdir + 'alignment.bam')

# def find_genes(alignment, ref, cov_threshold, fold_threshold, minbaseq, quiet=False):
#     # get length of each gene in ref
#     if not quiet: print("calculating gene lengths...")
#     gene_length = {}
#     total_length = 0
#     all_clusters = []
#     all_clusters_set = set()
#     for name, seq in pyfastx.Fasta(ref, build_index=False):
#         gene_length[name] = len(seq)
#         total_length += gene_length[name]
#         cluster = "__".join(name.split('__')[:2])
#         if cluster not in all_clusters_set:
#             all_clusters.append(cluster)
#         all_clusters_set.add(cluster)

#     # iterate over statistics, one record at a time
#     if not quiet: print("loading pileup stats...")
#     ref_support = Counter()
#     mybam = pysam.AlignmentFile(alignment)
#     for rec in tqdm(pysamstats.stat_variation(mybam, ref,
#         min_baseq=13), total=total_length):
#         if rec['matches'] >= fold_threshold:
#             ref_support[rec['chrom']] += 1

#     # collect genes with sufficient coverage
#     # [clusterUniqueIdentifier]__[clusterSymbol]__[alleleSymbol]__[alleleUniqueIdentifier]
#     if not quiet: print("filtering on coverage...")
#     found_clusters = defaultdict(list)
#     for gene in gene_length:
#         cov = ref_support[gene]/float(gene_length[gene])
#         print(ref_support[gene], gene_length[gene], cov)
#         if cov >= 0.5:
#             cluster = "__".join(gene.split('__')[:2])
#             found_clusters[cluster].append((cov, gene))

#     for cluster in found_clusters:
#         print(cluster)
#         found_clusters[cluster].sort(reverse=True)

#     return found_clusters, all_clusters

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args)

def mpile(seqs):
    if seqs is None: return([])
    thrbuf = mp.ThreadBuffer()
    return([(hit.ctg, hit.r_st-1, hit.r_en) for 
        hit in a.map(seqs[1], buf=thrbuf) if hit.NM<=3])

def generate_coverage(read1, read2, ref, pwid=0.98, ncpu=1, chunk_size=500000, quiet=FALSE):

    if not quiet: print("Building index and data structures...")

    seq_cov = {}
    for name, seq in pyfastx.Fasta(ref, build_index=False):
        seq_cov[name] = np.zeros(len(seq), dtype=int)

    nreads=0
    read_len
    for r in mp.fastx_read(read1):
        nreads+=1
        read_len = len(r[1])
    read_len /= nreads
    min_chain_score = int(pwid*read_len)

    a = mp.Aligner(ref, preset='sr', n_threads=ncpu, best_n=100000, min_chain_score=min_chain_score)  # load or build index
    if not a: raise Exception("ERROR: failed to load/build index")

    if not quiet: print("Aligning reads...")
    pool = ThreadPool(ncpu)
    for reads in tqdm(grouper(chain(
        mp.fastx_read(read1),
        mp.fastx_read(read2)), chunk_size), 
        total=int(1+2*nreads/chunk_size), disable=quiet):
        hits = pool.map(mpile, reads)
        for hit in chain.from_iterable(hits):
            if hit is None: continue
            seq_cov[hit[0]][hit[1]:hit[2]] += 1

    #close the pool and wait for the work to finish
    pool.close()
    pool.join()

    return(seq_cov)

def del_dups(seq):
    seen = set()
    pos = 0
    for item in seq:
        if item not in seen:
            seen.add(item)
            seq[pos] = item
            pos += 1
    del seq[pos:]
    return (seq)

def find_genes(coverage, mapping, cov_threshold, 
        prefix, outdir, fold_threshold=None, quiet=False):

    # get the best hit for each gene cluster (assumes srst2 format)
    best_hits = {}
    for gene in coverage:
        name = mapping[gene]
        cluster = tuple(name.split('__')[:2])
        if cluster in best_hits:
            if np.median(coverage[gene]) <= np.median(best_hits[cluster][1]):
                continue
        best_hits[cluster] = (name, coverage[gene])

    # get all clusters
    all_clusters = []
    for hit in mapping.values():
        all_clusters.append(tuple(hit.split('__')[:2]))
    all_clusters = del_dups(all_clusters)

    # output coverage information
    output_file = outdir + prefix + "_gene_coverage.csv"
    with open(output_file, 'w') as outfile:
        outfile.write("gene,position,depth\n")
        for cluster in best_hits:
            for pos,depth in enumerate(best_hits[cluster][1]):
                outfile.write(','.join([cluster, str(pos), str(depth)]) + '\n')

    # set fold threshold if not given
    fold_threshold = []
    for cluster in best_hits:
        fold_threshold.append(np.mean(best_hits[cluster][1][best_hits[cluster][1]>0]))
    fold_threshold = 0.1*np.mean(fold_threshold)

    # output cluster presence/absence
    output_file = outdir + prefix + "_gene_pa.csv"
    with open(output_file, 'w') as outfile:
        outfile.write("gene,found,coverage\n")
        for cluster in all_clusters:
            if cluster in best_hits:
                cov = np.sum(best_hits[cluster][1]>=fold_threshold)/float(len(best_hits[cluster][1]))
                if cov > cov_threshold:
                    outfile.write(cluster + ',1,' + str(cov) + '\n')
                else:
                    outfile.write(cluster + ',0,' + str(cov) + '\n')
            else:
                outfile.write(cluster + ',0,0\n')
    
    return


def get_options():
    import argparse

    description = 'Calls gene presence/absence using coverage and identity thresholds.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='call_genes_from_reads')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "--r1",
        dest="r1",
        required=True,
        help="Location of forward reads",
        type=str)
    io_opts.add_argument(
        "--r2",
        dest="r2",
        required=True,
        help="Location of reverse reads",
        type=str)
    io_opts.add_argument(
        "--db",
        dest="db",
        required=True,
        help="Fasta file of genes to be searched. Groups given in srst2 format.",
        type=str)

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of a new output directory",
                         type=str)

    filter_opts = parser.add_argument_group('Filter options')


    filter_opts.add_argument("--min_cov",
                         dest="min_cov",
                         help="minimum coverage to call a gene (default=0.9)",
                         type=float,
                         default=0.99)
    filter_opts.add_argument("--min_fold",
                         dest="min_fold",
                         help="read depth matching reference required to as 'covered' (default=2)",
                         type=float,
                         default=2)
    filter_opts.add_argument("--min_pwid",
                         dest="pwid",
                         help="minimum identity between read and reference (default=0.98)",
                         type=float,
                         default=0.98)


    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False)

    args = parser.parse_args()
    return (args)


def main():
    args = get_options()

    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # check files exist

    # clean database and remove sequences shorter than 300bp.
    temp_db = temp_dir + "temp_db.fasta"
    mapping = {}
    with open(temp_db, 'w') as outfile:
        index = 0
        for name, seq in pyfastx.Fasta(args.db, build_index=False):
            if len(seq)<=100: continue 
            outfile.write('>' + str(index) + '\n' + seq + '\n')
            mapping[str(index)] = name
            index += 1

    # align reads
    coverage = generate_coverage(
            read1=args.r1, 
            read2=args.r2,
            ref=temp_db,
            pwid=args.pwid, 
            ncpu=args.n_cpu, 
            chunk_size=500000, 
            quiet=args.quiet)

    # call genes
    found_genes, all_clusters = find_genes(coverage=coverage,
        mapping=mapping,
        cov_threshold=args.min_cov, 
        fold_threshold=args.min_fold,
        quiet=args.quiet)

    # write output
    prefix = os.path.basename(args.r1).split('.')[0].strip('_1')
    prefix += '_' + os.path.basename(args.db).split('.')[0]
    save_results(found_genes, all_clusters, prefix, args.output_dir)

    # clean up
    shutil.rmtree(temp_dir)
    
    return


if __name__ == '__main__':
    main()
