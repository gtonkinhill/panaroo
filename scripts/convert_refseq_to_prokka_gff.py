import sys, os
import argparse
import gffutils as gff
from io import StringIO
from Bio import SeqIO

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


def convert(gfffile, outputfile):

    #Split file and parse
    with open(gfffile, 'r') as infile:
        lines = infile.read()
    split = lines.split('##FASTA')

    if len(split) != 2:
        print("Problem reading GFF3 file: ", gff_file.name)
        raise RuntimeError("Error reading GFF3 input!")

    with StringIO(split[1]) as temp_fasta:
        sequences = list(SeqIO.parse(temp_fasta, 'fasta'))

    parsed_gff = gff.create_db(clean_gff_string(split[0]),
                               dbfn=":memory:",
                               force=True,
                               keep_order=False,
                               merge_strategy="create_unique",
                               from_string=True)

    with open(outputfile, 'w') as outfile:
        # write gff part
        outfile.write("##gff-version 3\n")
        for seq in sequences:
            outfile.write(" ".join(["##sequence-region", seq.id, "1", str(len(seq.seq))]) + "\n")
            
        prev_chrom=""
        prev_end=-1
        ids = set()
        for entry in parsed_gff.all_features(featuretype=()):
            # skip non CDS and overlapping regions
            if "CDS" not in entry.featuretype: continue
            if (entry.chrom==prev_chrom) and (entry.start<prev_end): continue
            # skip CDS that dont appear to be complete or have a premature stop codon
            
            premature_stop = False
            for sequence_index in range(len(sequences)):
                scaffold_id = sequences[sequence_index].id
                if scaffold_id == entry.seqid:
                    gene_sequence = sequences[sequence_index].seq[(entry.start -
                                                                1):entry.stop]
                    if (len(gene_sequence)%3 > 0) or (len(gene_sequence)<33): 
                        premature_stop = True
                        break
                    if entry.strand == "-":
                        gene_sequence = gene_sequence.reverse_complement()
                    if "*" in str(gene_sequence.translate())[:-1]:
                        premature_stop=True
                        break
            if premature_stop: continue
            
            c=1
            while entry.attributes['ID'][0] in ids:
                entry.attributes['ID'][0] += "." + str(c)
                c += 1
            ids.add(entry.attributes['ID'][0])
            prev_chrom = entry.chrom
            prev_end = entry.end
            print(entry, file=outfile)

        # write fasta part
        outfile.write("##FASTA\n")
        SeqIO.write(sequences, outfile, "fasta")

    return

def main():

    parser = argparse.ArgumentParser(description='Converts refseq GFF3 to prokka format.')

    parser.add_argument('-g', '--gff', dest='gff', type=str, required=True,
                       help='input gff file name')

    parser.add_argument('-o', '--out', dest='out', type=str, required=True,
                       help='output file name')

    args = parser.parse_args()

    convert(args.gff, args.out)

    return



if __name__ == '__main__':
    main()
