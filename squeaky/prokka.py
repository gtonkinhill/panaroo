#Takes .gff output from prokka and outputs combined gene/protien sequences from each isolate

import os

import gffutils as gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

#Clean other "##" starting lines from gff file, as it confuses parsers
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

#Get gene sequences from gff files
def get_gene_sequences(gff_list):
    sequence_list = []
    for gff_file in gff_list:
        #Get name and separate the prokka GFF into separate GFF and FASTA files
        name = gff_file.split('/')[-1].split('.')[0]
        temp_gff_name = 'temp_'+name+'.gff'
        temp_db_name = 'temp_'+name+".db"
        temp_fasta_name = 'temp_'+name+'.fasta'
        with open(gff_file) as temp_gff:
            lines = temp_gff.read()
            split = lines.split('##FASTA')
            with open(temp_fasta_name, 'w+') as temp_fasta:
                temp_fasta.write(split[1])
            with open(temp_gff_name, 'w+') as temp_gff:       
                temp_gff.write(clean_gff_string(split[0]))
        sequences = SeqIO.parse(temp_fasta_name, 'fasta')
        parsed_gff = gff.create_db(temp_gff_name, dbfn=temp_gff_name+'.db', force=True, keep_order=True)
        for entry in parsed_gff.all_features(featuretype=()):
            for sequence in sequences:
                if sequence.id == entry.seqid:
                    gene_sequence = sequence.seq[(entry.start-1):entry.stop]
            try:
                gene_name=entry.attributes["gene"]
                print(type(gene_name))
                print(gene_name)
            except KeyError:
                gene_name = "" 
            gene_record = SeqRecord(gene_sequence, id=entry.id, description=entry.attributes["product"], name=gene_name)
            sequence_list.append(gene_record)
    return sequence_list

#Translate sequences and return a second dic of protien sequences
def translate_sequences(gff_dic):
    protien_dic = {}
    for strain_id in gff_dic:
        sequences = None 

def output_sequence_dictionaries(protien_dictoinary, dna_dictionary):
    #Write the stuff to files combined_protien_CDS.fasta combined_DNA_CDS.fasta
    None

if __name__ == "__main__":
    #used for debugging purpopses
    import sys
    thing = get_gene_sequences([sys.argv[1]])
    print(len(thing))
