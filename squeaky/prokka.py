#Takes .gff output from prokka and outputs combined gene/protien sequences from each isolate

import os
from collections import OrderedDict
import gffutils as gff
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO


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


def get_gene_sequences(gff_file, file_number):
    #Get name and separate the prokka GFF into separate GFF and FASTA files
    sequence_dictionary = OrderedDict()

    #Split file and parse
    lines = gff_file.read()
    split = lines.split('##FASTA')
    with StringIO(split[1]) as temp_fasta:
        sequences = list(SeqIO.parse(temp_fasta, 'fasta'))

    parsed_gff = gff.create_db(
        clean_gff_string(split[0]),
        dbfn=":memory:",
        force=True,
        keep_order=True,
        from_string=True)

    #Get genes per scaffold
    scaffold_genes = {}
    for entry in parsed_gff.all_features(featuretype=()):
        if "CDS" not in entry.featuretype:
            continue
        scaffold_id = None
        for sequence_index in range(len(sequences)):
            scaffold_id = sequences[sequence_index].id
            if scaffold_id == entry.seqid:
                gene_sequence = sequences[sequence_index].seq[(
                    entry.start - 1):entry.stop]
                if entry.strand == "-":
                    gene_sequence = gene_sequence.reverse_complement()
                try:
                    gene_name = entry.attributes["gene"][0]
                except KeyError:
                    gene_name = " "
                gene_description = ";".join(entry.attributes["product"])
                gene_record = (entry.start,
                               SeqRecord(
                                   gene_sequence,
                                   id=entry.id,
                                   description=gene_description,
                                   name=gene_name,
                                   annotations={"scaffold":scaffold_id}))
                scaffold_genes[scaffold_id] = scaffold_genes.get(
                    scaffold_id, [])
                scaffold_genes[scaffold_id].append(gene_record)
    for scaffold in scaffold_genes:
        scaffold_genes[scaffold] = sorted(
            scaffold_genes[scaffold], key=lambda x: x[0])
    scaff_count = -1
    for scaffold in scaffold_genes:
        scaff_count += 1
        for gene_index in range(len(scaffold_genes[scaffold])):
            clustering_id = str(file_number) + '_' + str(
                scaff_count) + '_' + str(gene_index)
            sequence_dictionary[clustering_id] = scaffold_genes[scaffold][
                gene_index][1]

    return sequence_dictionary


#Translate sequences and return a second dic of protien sequences
def translate_sequences(sequence_dic):
    protein_list = []
    for strain_id in sequence_dic:
        sequence_record = sequence_dic[strain_id]
        protien_sequence = sequence_record.seq.translate()
        if protien_sequence[-1] == "*":
            protien_sequence = protien_sequence[0:-1]
        if "*" in protien_sequence:
            print(sequence_record)
            print(protien_sequence)
            raise ValueError("Premature stop codon in a gene!")
        protein_record = SeqRecord(
            protien_sequence, id=strain_id, description=strain_id)
        protein_list.append(protein_record)
    return protein_list


def output_files(dna_dictionary, prot_handle, dna_handle, csv_handle,
    gff_filename):
    #Simple output for protien list
    protien_list = translate_sequences(dna_dictionary)
    SeqIO.write(protien_list, prot_handle, 'fasta')
    #Correct DNA ids to CD-Hit acceptable ids, and output
    clustering_id_records = []
    for clusteringid in dna_dictionary:
        clean_record = SeqRecord(
            dna_dictionary[clusteringid].seq,
            id=clusteringid,
            description=clusteringid)
        clustering_id_records.append(clean_record)
    SeqIO.write(clustering_id_records, dna_handle, 'fasta')
    gff_name = os.path.splitext(os.path.basename(gff_filename.name))[0]

    #Combine everything to a csv and output it
    for protien in protien_list:
        clustering_id = protien.id
        relevant_seqrecord = dna_dictionary[clustering_id]
        out_list = [
            gff_name, relevant_seqrecord.annotations["scaffold"],
            clustering_id, relevant_seqrecord.id,
            str(protien.seq),
            str(relevant_seqrecord.seq),
            str(relevant_seqrecord.name),
            str(relevant_seqrecord.description)
        ]
        outline = ",".join(out_list)
        csv_handle.write(outline + '\n')
    return None


def process_prokka_input(gff_list, output_dir):
    try:
        protienHandle = open(output_dir + "combined_protein_CDS.fasta", 'w+')
        DNAhandle = open(output_dir + "combined_DNA_CDS.fasta", 'w+')
        csvHandle = open(output_dir + "gene_data.csv", 'w+')
        csvHandle.write(
            "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description\n"
        )
        for gff_no, gff in enumerate(gff_list):
            gene_sequences = get_gene_sequences(gff, gff_no)
            output_files(gene_sequences, protienHandle, DNAhandle, csvHandle, gff)
        protienHandle.close()
        DNAhandle.close()
        csvHandle.close()
        return True
    except:
        print("Error reading prokka input!")
        raise RuntimeError("Error reading prokka input!")


if __name__ == "__main__":
    #used for debugging purpopses
    import sys
    thing = process_prokka_input(sys.argv[1:])
