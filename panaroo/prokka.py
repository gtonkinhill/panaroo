#Takes .gff output from prokka and outputs combined gene/protien sequences from each isolate

import os
from collections import OrderedDict
import gffutils as gff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import generic_by_id
from Bio.SeqRecord import SeqRecord
from io import StringIO
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
from .biocode_convert import convert_gbk_gff3

bact_translation_table = np.array([[[b'K', b'N', b'K', b'N', b'X'],
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

def get_trans_table(table):
    # swap to different codon table
    translation_table = bact_translation_table.copy()
    tb = generic_by_id[table]
    if table!=11:
        if table not in generic_by_id:
            raise RuntimeError("Invalid codon table! Must be available" +
                " as a generic table in BioPython")
        for codon in tb.forward_table:
            if 'U' in codon: continue
            ind = reduce_array[np.array(bytearray(codon.encode()), dtype=np.int8)]
            translation_table[ind[0], ind[1], ind[2]] = tb.forward_table[codon].encode('utf-8')
        for codon in tb.stop_codons:
            if 'U' in codon: continue
            ind = reduce_array[np.array(bytearray(codon.encode()), dtype=np.int8)]
            translation_table[ind[0], ind[1], ind[2]] = b'*'

    return([translation_table, set(tb.start_codons)])


def translate(seq, translation_table):
    indices = reduce_array[np.array(bytearray(seq.encode()), dtype=np.int8)]
    pseq = translation_table[0][
        indices[np.arange(0, len(seq), 3)], indices[np.arange(1, len(seq), 3)],
        indices[np.arange(2, len(seq), 3)]].tostring().decode('ascii')
    # Check for a different start codon.
    if seq[0:3] in translation_table[1]:
        return ('M' + pseq[1:])
    return(pseq)


def create_temp_gff3(gff_file, fasta_file, temp_dir):

    # create directory if it isn't present already
    if not os.path.exists(temp_dir + "temp_gffs"):
        os.mkdir(temp_dir + "temp_gffs")
    
    prefix = os.path.splitext(os.path.basename(gff_file))[0]
    ext = os.path.splitext(gff_file)[1]

    if fasta_file is None:
        try:
            convert_gbk_gff3(gff_file, temp_dir + "temp_gffs/" + prefix + '.gff', True)
        except:
            print("Error reading Genbank input! These must compliant with Genbank/ENA/DDJB." +
                  " This can be forced in Prokka by specifying the --compliance parameter.")
            raise RuntimeError(
            "Error reading Genbank input: {0}\n"
            "These must compliant with Genbank/ENA/DDJB. This can be forced in Prokka"
            " by specifying the --compliance parameter."
            .format(gff_file)
            )
    else:
        if ext not in ['.gff', '.gff3']:
            raise RuntimeError(f"Invalid file extension! ({ext})")
        fasta_ext = os.path.splitext(fasta_file)[1]
        if fasta_ext not in [".fasta", ".fa", ".fas", ".fna"]:
            raise RuntimeError(f"Invalid file extension! ({fasta_ext})")

        # merge files into temporary gff3
        with open(temp_dir + "temp_gffs/" + prefix + '.gff', 'w') as outfile:
            with open(gff_file, 'r') as infile:
                gff_string = infile.read().strip()
                if '\naccn' in gff_string: # deal with PATRIC input format
                    gff_string = gff_string.replace('accn|', '')
                outfile.write(gff_string)
                outfile.write('\n##FASTA\n')
            with open(fasta_file, 'r') as infile:
                outfile.write(infile.read().strip())

    return(temp_dir + "temp_gffs/" + prefix + '.gff')

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


def get_gene_sequences(gff_file_name, file_number, filter_seqs, table):
    #Get name and separate the prokka GFF into separate GFF and FASTA files
    if ',' in gff_file_name:
        print("Problem reading GFF3 file: ", gff_file_name)
        raise RuntimeError("Error reading prokka input!")

    gff_file = open(gff_file_name, 'r')
    sequence_dictionary = OrderedDict()

    #Split file and parse
    lines = gff_file.read().replace(',', '')
    split = lines.split('##FASTA')

    if len(split) != 2:
        print("Problem reading GFF3 file: ", gff_file.name)
        raise RuntimeError("Error reading prokka input!")

    with StringIO(split[1]) as temp_fasta:
        sequences = list(SeqIO.parse(temp_fasta, 'fasta'))

    parsed_gff = gff.create_db(clean_gff_string(split[0]),
                               dbfn=":memory:",
                               force=True,
                               keep_order=True,
                               from_string=True,
                               merge_strategy="create_unique")

    #Get genes per scaffold
    scaffold_genes = {}
    for entry in parsed_gff.all_features(featuretype=()):
        if "CDS" not in entry.featuretype:
            continue

        scaffold_id = None
        gene_sequence = None
        for sequence_index in range(len(sequences)):
            scaffold_id = sequences[sequence_index].id
            if scaffold_id == entry.seqid:
                gene_sequence = sequences[sequence_index].seq[(entry.start -
                                                               1):entry.stop]
                if entry.strand == "-":
                    gene_sequence = gene_sequence.reverse_complement()
                try:
                    gene_name = entry.attributes["gene"][0]
                except KeyError:
                    gene_name = ""
                if gene_name == "":
                    try:
                        gene_name = entry.attributes["name"][0]
                    except KeyError:
                        gene_name = ""

                try:
                    gene_description = ";".join(entry.attributes["product"])
                    gene_description = gene_description.replace(",", "")
                except KeyError:
                    gene_description = ""

                #clean entries if requested
                if entry.frame != '0':
                    print('Invalid gene! Panaroo currently does not support frame shifts.')
                    if filter_seqs: continue
                    else: raise ValueError("Invalid gene sequence!")

                if ((len(gene_sequence) % 3 > 0) or
                    (len(gene_sequence) < 34)) or ("*" in translate(str(gene_sequence), table)[:-1]):
                    print('invalid gene! file - id: ', gff_file_name, ' - ',
                          entry.id)
                    print('Length:', len(gene_sequence), ', Has stop:', ("*" in str(
                        gene_sequence.translate())[:-1]))
                    if filter_seqs: continue
                    else: raise ValueError("Invalid gene sequence!")

                gene_record = (entry.start,
                               SeqRecord(gene_sequence,
                                         id=entry.id,
                                         description=gene_description,
                                         name=gene_name,
                                         annotations={"scaffold":
                                                      scaffold_id}))
                scaffold_genes[scaffold_id] = scaffold_genes.get(
                    scaffold_id, [])
                scaffold_genes[scaffold_id].append(gene_record)
        if gene_sequence is None:
            print('Sequence ID not found in Fasta!', entry.seqid)
            if filter_seqs: continue
            else: raise ValueError("Invalid gene sequence!")

    if len(scaffold_genes) == 0:
        print("No valid sequences found in GFF!", gff_file_name)
        raise ValueError("Invalid GFF!")

    for scaffold in scaffold_genes:
        scaffold_genes[scaffold] = sorted(scaffold_genes[scaffold],
                                          key=lambda x: x[0])
    scaff_count = -1
    for scaffold in scaffold_genes:
        scaff_count += 1
        for gene_index in range(len(scaffold_genes[scaffold])):
            clustering_id = str(file_number) + '_' + str(
                scaff_count) + '_' + str(gene_index)
            sequence_dictionary[clustering_id] = scaffold_genes[scaffold][
                gene_index][1]

    gff_file.close()

    return sequence_dictionary, translate_sequences(sequence_dictionary, table)


#Translate sequences and return a second dic of protien sequences
def translate_sequences(sequence_dic, table):
    protein_list = []
    for strain_id in sequence_dic:
        sequence_record = sequence_dic[strain_id]
        if (len(sequence_record.seq) % 3) != 0:
            raise ValueError(
                "Coding sequence not divisible by 3, is it complete?!")
        protien_sequence = translate(str(sequence_record.seq), table)
        if protien_sequence[-1] == "*":
            protien_sequence = protien_sequence[0:-1]
        if "*" in protien_sequence:
            print(sequence_record)
            print(protien_sequence)
            raise ValueError("Premature stop codon in a gene!")
        protein_record = SeqRecord(Seq(protien_sequence),
                                   id=strain_id,
                                   description=strain_id)
        protein_list.append(protein_record)
    return protein_list


def output_files(dna_dictionary, protien_list, prot_handle, dna_handle,
                 csv_handle, gff_filename):
    #Simple output for protien list
    SeqIO.write(protien_list, prot_handle, 'fasta')
    #Correct DNA ids to CD-Hit acceptable ids, and output
    clustering_id_records = []
    for clusteringid in dna_dictionary:
        clean_record = SeqRecord(dna_dictionary[clusteringid].seq,
                                 id=clusteringid,
                                 description=clusteringid)
        clustering_id_records.append(clean_record)
    SeqIO.write(clustering_id_records, dna_handle, 'fasta')
    gff_name = os.path.splitext(os.path.basename(gff_filename))[0]

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


def process_prokka_input(gff_list, output_dir, filter_seqs, quiet, n_cpu, table):
    trans_table = get_trans_table(table)
    try:
        protienHandle = open(output_dir + "combined_protein_CDS.fasta", 'w+')
        DNAhandle = open(output_dir + "combined_DNA_CDS.fasta", 'w+')
        csvHandle = open(output_dir + "gene_data.csv", 'w+')
        csvHandle.write(
            "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description\n"
        )
        job_list = list(enumerate(gff_list))
        job_list = [
            job_list[i:i + n_cpu] for i in range(0, len(job_list), n_cpu)
        ]
        for job in tqdm(job_list, disable=quiet):
            gene_sequence_list = Parallel(n_jobs=n_cpu)(
                delayed(get_gene_sequences)(gff, gff_no, filter_seqs, trans_table)
                for gff_no, gff in job)
            for i, gene_seq in enumerate(gene_sequence_list):
                output_files(gene_seq[0], gene_seq[1], protienHandle,
                             DNAhandle, csvHandle, job[i][1])
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
    thing = process_prokka_input([open(f, 'r') for f in sys.argv[1:]], "./", 2)
