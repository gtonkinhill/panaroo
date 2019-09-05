import sys, os
import argparse
from collections import OrderedDict, defaultdict
import gffutils as gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import numpy as np
from random import sample
from dendropy.simulate import treesim
from dendropy.model import reconcile
from dendropy import TaxonNamespace
import copy
import math

codons = [
    'ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 'ACT', 'AAC', 'AAT',
    'AAA', 'AAG', 'AGC', 'AGT', 'AGA', 'AGG', 'CTA', 'CTC', 'CTG', 'CTT',
    'CCA', 'CCC', 'CCG', 'CCT', 'CAC', 'CAT', 'CAA', 'CAG', 'CGA', 'CGC',
    'CGG', 'CGT', 'GTA', 'GTC', 'GTG', 'GTT', 'GCA', 'GCC', 'GCG', 'GCT',
    'GAC', 'GAT', 'GAA', 'GAG', 'GGA', 'GGC', 'GGG', 'GGT', 'TCA', 'TCC',
    'TCG', 'TCT', 'TTC', 'TTT', 'TTA', 'TTG', 'TAC', 'TAT', 'TGC', 'TGT', 'TGG'
]
codons = [Seq(c) for c in codons]

translation_table = np.array([[[b'K', b'N', b'K', b'N', b'X'],
                               [b'T', b'T', b'T', b'T', b'T'],
                               [b'R', b'S', b'R', b'S', b'X'],
                               [b'I', b'I', b'M', b'I', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'Q', b'H', b'Q', b'H', b'X'],
                               [b'P', b'P', b'P', b'P', b'P'],
                               [b'R', b'R', b'R', b'R', b'R'],
                               [b'L', b'L', b'L', b'L', b'L'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'E', b'D', b'E', b'D', b'X'],
                               [b'A', b'A', b'A', b'A', b'A'],
                               [b'G', b'G', b'G', b'G', b'G'],
                               [b'V', b'V', b'V', b'V', b'V'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'*', b'Y', b'*', b'Y', b'X'],
                               [b'S', b'S', b'S', b'S', b'S'],
                               [b'*', b'C', b'W', b'C', b'X'],
                               [b'L', b'F', b'L', b'F', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']],
                              [[b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X'],
                               [b'X', b'X', b'X', b'X', b'X']]])

reduce_array = np.full(200, 4)
reduce_array[[65, 97]] = 0
reduce_array[[67, 99]] = 1
reduce_array[[71, 103]] = 2
reduce_array[[84, 116]] = 3


def translate(seq):

    indices = reduce_array[np.fromstring(seq, dtype=np.int8)]

    return translation_table[indices[np.arange(0, len(
        seq), 3)], indices[np.arange(1, len(seq), 3)], indices[np.arange(
            2, len(seq), 3)]].tostring().decode('ascii')


def get_codon(index, strand="+"):
    codon = codons[index]
    if strand == "-":
        codon = codon.reverse_complement()
    return np.array(list(str(codon)))


def clean_gff_string(gff_string):
    splitlines = gff_string.splitlines()
    lines_to_delete = []
    for index in range(len(splitlines)):
        if '##sequence-region' in splitlines[index]:
            lines_to_delete.append(index)
    for index in sorted(lines_to_delete, reverse=True):
        del splitlines[index]
    cleaned_gff = "\n".join(splitlines)
    return cleaned_gff


def simulate_img_with_mutation(in_tree,
                               gain_rate,
                               loss_rate,
                               mutation_rate,
                               ngenes=100,
                               min_ncore=10):
    # simulate accessory p/a using infintely many genes model
    n_additions = 0
    for node in in_tree.preorder_node_iter():
        node.acc_genes = []
        if node.parent_node is not None:
            # simulate loss of genes from previous node
            node.acc_genes = [
                g for g in node.acc_genes
                if np.random.poisson(lam=node.edge.length * loss_rate / 2.0,
                                     size=1) > 0
            ]
            # simulate new genes with lengths sampled uniformly.
            n_new = np.random.poisson(lam=node.edge.length * gain_rate / 2.0,
                                      size=1)[0]
            lengths = np.random.uniform(low=0.0,
                                        high=node.edge.length,
                                        size=n_new)
            for l in lengths:
                # simulate loss using this length
                if np.random.poisson(lam=l * loss_rate / 2.0, size=1)[0] > 0:
                    n_new -= 1
            # add new genes to node
            node.acc_genes = node.parent_node.acc_genes + list(
                range(n_additions, n_additions + n_new))
            n_additions += n_new

    print("accessory size: ", n_additions)

    # Now add core
    ncore = ngenes - n_additions
    if ncore < min_ncore:
        ncore = min_ncore

    core_genes = list(range(n_additions, n_additions + ncore))
    for node in in_tree.preorder_node_iter():
        node.acc_genes += core_genes

    # Now add mutations
    n_condons = len(codons)
    for node in in_tree.preorder_node_iter():
        node.gene_mutations = defaultdict(list)
        if node.parent_node is not None:
            # copy mutations from parent
            for g in node.acc_genes:
                if g in node.parent_node.gene_mutations:
                    node.gene_mutations[g] = node.parent_node.gene_mutations[
                        g].copy()
            # add mutations
            for g in node.acc_genes:
                n_new = np.random.poisson(lam=node.edge.length *
                                          mutation_rate / 2.0,
                                          size=1)[0]
                locations = list(np.random.uniform(low=0.0, high=1,
                                                   size=n_new))
                mutations = [(sample(range(0, n_condons), 1)[0], l)
                             for l in locations]
                node.gene_mutations[g] += mutations

    return in_tree


def simulate_pangenome(ngenes, nisolates, effective_pop_size, gain_rate,
                       loss_rate, mutation_rate):

    # simulate a phylogeny using the coalscent
    sim_tree = treesim.pure_kingman_tree(taxon_namespace=TaxonNamespace(
        [str(i) for i in range(1, 1 + nisolates)]),
                                         pop_size=effective_pop_size)

    # simulate gene p/a and mutation
    sim_tree = simulate_img_with_mutation(sim_tree,
                                          gain_rate=gain_rate,
                                          loss_rate=loss_rate,
                                          mutation_rate=mutation_rate,
                                          ngenes=ngenes)

    # get genes and mutations for each isolate
    gene_mutations = []
    for leaf in sim_tree.leaf_node_iter():
        gene_mutations.append([[g, leaf.gene_mutations[g]]
                               for g in leaf.acc_genes])

    return (gene_mutations)


def add_diversity(gfffile, nisolates, effective_pop_size, gain_rate, loss_rate,
                  mutation_rate, n_sim_genes, prefix):

    with open(gfffile, 'r') as infile:
        lines = infile.read()

    split = lines.split('##FASTA')
    if len(split) != 2:
        print("Problem reading GFF3 file: ", gfffile)
        raise RuntimeError("Error reading GFF3 input!")

    with StringIO(split[1]) as temp_fasta:
        sequences = list(SeqIO.parse(temp_fasta, 'fasta'))
    seq_dict = OrderedDict()
    for seq in sequences:
        seq_dict[seq.id] = np.array(list(str(seq.seq)))

    parsed_gff = gff.create_db(clean_gff_string(split[0]),
                               dbfn=":memory:",
                               force=True,
                               keep_order=False,
                               merge_strategy="create_unique",
                               sort_attribute_values=True,
                               from_string=True)

    #Get gene entries to modify
    all_gene_locations = []
    gene_locations = []
    prev_end = -1
    gene_seqs = []

    for entry in parsed_gff.all_features(featuretype=()):
        if "CDS" not in entry.featuretype: continue

        left = entry.start - 1
        right = entry.stop
        gene_sequence = Seq(''.join(seq_dict[entry.seqid][left:right]))
        if entry.strand == "-":
            gene_sequence = gene_sequence.reverse_complement()
        gene_sequence = gene_sequence.translate()
        gene_seqs.append(SeqRecord(gene_sequence, id=entry.id, description=""))

        all_gene_locations.append(entry)
        if entry.start < prev_end:
            prev_end = entry.end
            gene_locations = gene_locations[0:-1]
            continue
        prev_end = entry.end
        gene_locations.append(entry)

    # sub-sample genes so that some are conserved
    gene_locations = sample(gene_locations, n_sim_genes)

    # simulate presence/absence matrix and gene mutations (only swap codons)
    pan_sim = simulate_pangenome(ngenes=len(gene_locations),
                                 nisolates=nisolates,
                                 effective_pop_size=effective_pop_size,
                                 gain_rate=gain_rate,
                                 loss_rate=loss_rate,
                                 mutation_rate=mutation_rate)

    #Modify each gene
    for i, pan in enumerate(pan_sim):
        temp_seq_dict = copy.deepcopy(seq_dict)
        included_genes = set()
        n_mutations = 0
        for gene in pan:
            entry = gene_locations[gene[0]]
            included_genes.add(gene[0])

            left = entry.start - 1
            right = entry.stop
            if right < left: raise RuntimeError("Error issue with left/right!")

            start_sites = list(range(left, right, 3))[1:-1]

            n_mutations += len(gene[1])
            # swap codons at chosen start sites
            for mutation in gene[1]:
                # find start site of codon swap
                start = start_sites[math.floor(mutation[1] * len(start_sites))]
                cod = get_codon(index=mutation[0], strand=entry.strand)
                if (start < left) or ((start + 3) > (right)):
                    raise RuntimeError("Error issue with start!")
                temp_seq_dict[entry.seqid][start:(start + 3)] = cod

        # remove genes not in the accessory
        deleted_genes = 0
        d_index = defaultdict(lambda: np.array([]))
        for g, entry in enumerate(gene_locations):
            left = entry.start - 1
            right = entry.stop
            if right < left: raise RuntimeError("Error issue with left/right!")
            if g not in included_genes:
                deleted_genes += 1
                d_index[entry.seqid] = np.append(d_index[entry.seqid],
                                                 np.arange(left, right))

            gene_sequence = Seq(''.join(
                temp_seq_dict[entry.seqid][left:right]))
            if entry.strand == "-":
                gene_sequence = gene_sequence.reverse_complement()
            gene_sequence = gene_sequence.translate()
            gene_seqs.append(
                SeqRecord(gene_sequence, id=entry.id, description=""))

        for entryid in d_index:
            temp_seq_dict[entryid] = np.delete(temp_seq_dict[entry.seqid],
                                               d_index[entryid])

        print("mutations in genome: ", n_mutations)
        print("genes deleted: ", deleted_genes)

        # write out sequences
        out_name = prefix + "_iso_" + str(i) + ".fasta"
        outfile = open(out_name, 'w')
        sequences = [
            SeqRecord(Seq(''.join(temp_seq_dict[s])), id=s, description="")
            for s in temp_seq_dict
        ]
        SeqIO.write(sequences, outfile, 'fasta')
        # close file
        outfile.close()

    # write out database for prokka
    prokka_db_name = prefix + "_prokka_DB.fasta"
    with open(prokka_db_name, 'w') as dboutfile:
        SeqIO.write(gene_seqs, dboutfile, 'fasta')

    # write presence/absence file
    pa_by_iso = []
    for i, pan in enumerate(pan_sim):
        pa = set()
        for gene in pan:
            pa.add(gene[0])
        pa_by_iso.append(pa)

    out_name = prefix + "_presence_absence.csv"

    seen = set()
    with open(out_name, 'w') as outfile:
        outfile.write("\t".join(
            ["Gene"] + ["iso" + str(i)
                        for i in range(1, nisolates + 1)]) + "\n")
        for g, entry in enumerate(gene_locations):
            seen.add(entry.id)
            outfile.write("\t".join(
                [entry.id] +
                ["1" if g in pa_by_iso[i] else "0"
                 for i in range(nisolates)]) + "\n")

        for g, entry in enumerate(all_gene_locations):
            if entry.id in seen: continue
            outfile.write("\t".join([entry.id] +
                                    ["1" for i in range(nisolates)]) + "\n")

    return


def main():

    parser = argparse.ArgumentParser(description=(
        'Simulates a pangenome using the infinitely many genes ' +
        'model and adds mutational variation to genes. Takes a gff3 file as input.'
    ))

    parser.add_argument('-g',
                        '--gff',
                        dest='gff',
                        type=str,
                        required=True,
                        help='input gff file name')

    parser.add_argument('--nisolates',
                        dest='nisolates',
                        type=int,
                        required=True,
                        help='number of genomes to simulate')

    parser.add_argument('--mutation_rate',
                        dest='mutation_rate',
                        type=float,
                        required=True,
                        help='mutation rate of genes')

    parser.add_argument('--gain_rate',
                        dest='gain_rate',
                        type=float,
                        required=True,
                        help='gain rate of accessory genes')

    parser.add_argument('--loss_rate',
                        dest='loss_rate',
                        type=float,
                        required=True,
                        help='loss rate of accessory genes')

    parser.add_argument('--pop_size',
                        dest='pop_size',
                        type=float,
                        required=True,
                        help='effective population size')

    parser.add_argument(
        '--n_sim_genes',
        dest='n_sim_genes',
        type=int,
        required=True,
        help=('max number of genes that may be ' +
              'affected by the simulation. The rest' + ' will be left as is.'))

    parser.add_argument('-o',
                        '--out',
                        dest='output_dir',
                        type=str,
                        required=True,
                        help='output directory')

    args = parser.parse_args()

    args.pop_size = math.floor(args.pop_size)
    args.output_dir = os.path.join(args.output_dir, "")

    prefix = (args.output_dir + "pan_sim_gr_" + str(args.gain_rate) + "_lr_" +
              str(args.loss_rate) + "_mu_" + str(args.mutation_rate))

    # adjust rates for popsize
    args.gain_rate = 2.0 * args.gain_rate * args.pop_size
    args.loss_rate = 2.0 * args.loss_rate * args.pop_size
    args.mutation_rate = 2.0 * args.mutation_rate * args.pop_size

    add_diversity(gfffile=args.gff,
                  nisolates=args.nisolates,
                  effective_pop_size=args.pop_size,
                  gain_rate=args.gain_rate,
                  loss_rate=args.loss_rate,
                  mutation_rate=args.mutation_rate,
                  n_sim_genes=args.n_sim_genes,
                  prefix=prefix)

    return


if __name__ == '__main__':
    main()
