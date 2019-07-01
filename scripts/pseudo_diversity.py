import sys, os
import argparse
from collections import OrderedDict
import gffutils as gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import numpy as np
from random import sample 

codons = ['ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 
    'ACT', 'AAC', 'AAT', 'AAA', 'AAG', 'AGC', 'AGT', 'AGA', 
    'AGG', 'CTA', 'CTC', 'CTG', 'CTT', 'CCA', 'CCC', 'CCG', 
    'CCT', 'CAC', 'CAT', 'CAA', 'CAG', 'CGA', 'CGC', 'CGG', 
    'CGT', 'GTA', 'GTC', 'GTG', 'GTT', 'GCA', 'GCC', 'GCG', 
    'GCT', 'GAC', 'GAT', 'GAA', 'GAG', 'GGA', 'GGC', 'GGG', 
    'GGT', 'TCA', 'TCC', 'TCG', 'TCT', 'TTC', 'TTT', 'TTA', 
    'TTG', 'TAC', 'TAT', 'TGC', 'TGT', 'TGG']
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

def random_codon(strand = "+"):
    # codon = sample(codons, 1)[0]
    codon = codons[1]
    if strand=="-":
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

def add_diversity(gfffile, outputfile, ngenes, diversity, sampled_entries=None):

    outfile = open(outputfile, 'w')

    with open(gfffile, 'r') as infile:
        lines = infile.read()

    split = lines.split('##FASTA')
    if len(split) != 2:
        print("Problem reading GFF3 file: ", gfffile)
        raise RuntimeError("Error reading GFF3 input!")

    outfile.write(split[0])
    outfile.write('##FASTA\n')
    
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
    if sampled_entries is None:
        sampled_entries = []
        prev_end = -1
        for entry in parsed_gff.all_features(featuretype=()):
            if "CDS" not in entry.featuretype: continue
            if entry.start < prev_end:
                prev_end = entry.end
                sampled_entries = sampled_entries[0:-1]
                continue
            prev_end = entry.end
            sampled_entries.append(entry)
        index = np.random.choice(len(sampled_entries), size=ngenes, replace=False)
        sampled_entries = [sampled_entries[i] for i in index]

    #Modify each gene
    for entry in sampled_entries:
        left = entry.start - 1
        right = entry.stop
        if right<left: raise RuntimeError("Error issue with left/right!")

        condon_start_sites = np.random.choice(range(left, right, 3), 
            size=int(np.ceil((1.0-diversity)*(right-left)/3.0)), 
            replace=True)
        # swap codons at chosen start sites
        for start in condon_start_sites:
            cod = random_codon(strand = entry.strand)
            if (start<left) or ((start+3)>(right)):
                raise RuntimeError("Error issue with start!")
            seq_dict[entry.seqid][start:(start+3)] = cod
        # if entry.strand=="-":
        #     prot = str(translate(str(Seq(''.join(seq_dict[entry.seqid][left:right])).reverse_complement())))
        # else:
        #     prot = str(translate(str(Seq(''.join(seq_dict[entry.seqid][left:right])))))
        # print(prot)
        # if "*" in prot[0:-1]:
        #     raise RuntimeError("Error stop where there shouldnt be!")

    # write out sequences
    sequences = [SeqRecord(Seq(''.join(seq_dict[s])),
                                id=s,
                                description="") for s in seq_dict]
    SeqIO.write(sequences, outfile, 'fasta')

    # close file
    outfile.close()

    return sampled_entries


def main():

    parser = argparse.ArgumentParser(
        description='Adds artificial variation to genes in a gff.')

    parser.add_argument('-g',
                        '--gff',
                        dest='gff',
                        type=str,
                        required=True,
                        help='input gff file name')

    parser.add_argument('--ngenes',
                        dest='ngenes',
                        type=int,
                        required=True,
                        help='number of genes to add diversity to')

    parser.add_argument('--nreps',
                        dest='nreps',
                        type=int,
                        required=True,
                        help='number of replications')

    parser.add_argument('-d',
                        '--diversity',
                        dest='diversity',
                        type=float,
                        required=True,
                        help='induced percentage diversity in genes')

    parser.add_argument('-o',
                        '--out',
                        dest='output_dir',
                        type=str,
                        required=True,
                        help='output directory')

    args = parser.parse_args()

    args.output_dir = os.path.join(args.output_dir, "")

    prefix = os.path.splitext(os.path.basename(args.gff))[0]
    sampled_entries = None
    for i in range(args.nreps):
        out_file_name = (args.output_dir + prefix + "div_" +
                        str(args.diversity) + "_rep_" + str(i) + ".gff")

        add_diversity(args.gff, out_file_name, args.ngenes, args.diversity,
            sampled_entries)

    return


if __name__ == '__main__':
    main()
