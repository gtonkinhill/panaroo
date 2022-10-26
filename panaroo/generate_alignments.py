import os
import subprocess
import sys
import re

import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from Bio.Align.Applications import PrankCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
import Bio.Application
from Bio import BiopythonExperimentalWarning
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import codonalign

unambiguous_degenerate_codons = {"ACN":"T", "TCN":"S", "CTN":"L", "CCN":"P",
                                 "CGN":"R", "GTN":"V", "GCN":"A", "GGN":"G"}


def check_aligner_install(aligner):
    """Checks for the presence of the specified aligned in $PATH

    Args:
        check_aligner_install(str)
            str = specified aligner

    Returns:
        presence (bool)
            True/False aligner present
    """
    if aligner == "clustal":
        command = "clustalo --help"
    elif aligner == "prank":
        command = "prank -help"
    elif aligner == "mafft":
        command = "mafft --help"
    else:
        sys.stderr.write("Incorrect aligner specification\n")
        sys.exit()
    p = str(
        subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
        )
    )
    present = False

    if aligner == "clustal":
        find_ver = re.search(r"Clustal Omega - \d+\.\d+\.\d+", p)
    elif aligner == "prank":
        find_ver = re.search(r"prank v\.\d+\.", p)
    elif aligner == "mafft":
        find_ver = re.search(r"MAFFT v\d+\.\d+", p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need specified aligner to be installed " + "\n")
        sys.exit(1)

    return present


def output_sequence(node, isolate_list, temp_directory, outdir):
    # Get the name of the sequences for the gene of interest
    sequence_ids = node["seqIDs"]
    output_sequences = []
    #Counter for the number of sequences for downstream check of >1
    isolate_no = 0
    # Look for gene sequences among all genes (from disk)
    for seq in SeqIO.parse(outdir + "combined_DNA_CDS.fasta", "fasta"):
        isolate_num = int(seq.id.split("_")[0])
        isolate_name = isolate_list[isolate_num].replace(";", "") + ";" + seq.id
        if seq.id in sequence_ids:
            output_sequences.append(SeqRecord(seq.seq, id=isolate_name, description=""))
            isolate_no += 1
    # Put gene of interest sequences in a generator, with corrected isolate names
    output_sequences = (x for x in output_sequences)
    # set filename to gene name, if more than one sequence to be aliged
    if isolate_no > 1:
        outname = temp_directory + node["name"] + ".fasta"
    else:
        # If only one sequence, output it to aliged directory and break
        outname = outdir + "/aligned_gene_sequences/" + node["name"] + ".fasta"
        SeqIO.write(output_sequences, outname, "fasta")
        return None
    # check to see if filename is too long
    if len(outname) >= 248:
        outname = outname[:248] + ".fasta"
    # Write them to disk
    SeqIO.write(output_sequences, outname, "fasta")
    return outname

def output_dna_and_protein(node, isolate_list, temp_directory, outdir, 
                           all_proteins, all_dna):
    #Get the name of the sequences for the gene of interest

    sequence_ids = node["seqIDs"]
    #Make sure later steps don't iterate through a string
    if type(sequence_ids) == str:
        sequence_ids = [sequence_ids]
    output_dna = []
    output_protein = []
    #Counter for the number of sequences to avoid aliging single sequences
    isolate_no = 0


    
    for seq_id in sequence_ids:
        #Get isolate names
        isolate_num = int(seq_id.split('_')[0])
        isolate_name = isolate_list[isolate_num].replace(";",
                                                         "") + ";" + seq_id
        #prot_stop_codon = all_proteins[seq_id].seq.find("*")
            
            
        #if prot_stop_codon > -1:
        #    output_protein.append(
        #        SeqRecord(all_proteins[seq_id].seq[:prot_stop_codon], 
        #                  id=isolate_name, description=""))
        #    output_dna.append(
        #        SeqRecord(all_dna[seq_id].seq[:(prot_stop_codon*3)], 
        #                  id=isolate_name, description=""))
        #else:
        output_dna.append(
            SeqRecord(all_dna[seq_id].seq, 
                      id=isolate_name, description=""))
        output_protein.append(
            SeqRecord(all_proteins[seq_id].seq, 
                      id=isolate_name, description=""))
        isolate_no += 1
    

    #only output genes with more than one isolate in them
    if isolate_no > 1:
        #set filename to gene name
        prot_outname = temp_directory + node["name"] + ".fasta"
        if len(prot_outname) > 248:
            prot_outname = prot_outname[:248] + '.fasta'
        dna_outname = outdir + "unaligned_dna_sequences/" + node["name"] + ".fasta"
        if len(dna_outname) > 248:
            dna_outname = dna_outname[:248] + '.fasta'
        #Write them to disk time
        SeqIO.write(output_protein, prot_outname, 'fasta')
        SeqIO.write(output_dna, dna_outname, 'fasta')
        output_files = (prot_outname, dna_outname)
        
    else:
        output_singleton = output_dna
        #set filename, write
        singleton_outname = outdir + "aligned_gene_sequences/" + node["name"] +".aln.fas"
        SeqIO.write(output_singleton, singleton_outname, 'fasta')
        output_files = (None, None)
        
    return output_files


def get_alignment_commands(fastafile_name, outdir, aligner, threads):
    geneName = fastafile_name.split("/")[-1].split(".")[0]
    if aligner == "prank":
        command = PrankCommandline(d=fastafile_name, 
            o=outdir + "aligned_gene_sequences/" + geneName, 
            f=8, codon=False
        )
    elif  aligner == "mafft":
        command = MafftCommandline(
            input=fastafile_name, auto=True, nuc=True, adjustdirection=True, thread=1
        )
    elif aligner == "clustal":
        command = ClustalOmegaCommandline(
            infile=fastafile_name,
            outfile=outdir + "aligned_gene_sequences/" + geneName + ".aln.fas",
            seqtype="DNA",  threads=1
        )
   
    return (command, fastafile_name)

def get_protein_commands(fastafile_name, outdir, aligner, threads):
    if fastafile_name != None:
        geneName = fastafile_name.split('/')[-1].split('.')[0]
    else:
        return (None, None)
    if aligner == "prank":
        command = PrankCommandline(d=fastafile_name,
                                   o=geneName,
                                   f=8,
                                   protein=True)
    elif (threads > 3):
        if aligner == "mafft":
            command = MafftCommandline(input=fastafile_name,
                                       auto=True,
                                       amino=True)
        elif aligner == "clustal":
            command = ClustalOmegaCommandline(
                infile=fastafile_name,
                outfile=outdir + "aligned_protein_sequences/" + geneName +
                ".aln.fas",
                seqtype="Protein")
    elif (threads <= 3):
        if aligner == "mafft":
            command = MafftCommandline(input=fastafile_name,
                                       auto=True,
                                       thread=threads,
                                       amino=True)
        elif aligner == "clustal":
            command = ClustalOmegaCommandline(
                infile=fastafile_name,
                outfile=outdir + "aligned_protein_sequences/" + geneName +
                ".aln.fas",
                seqtype="Protein",
                threads=threads)
    return (command, fastafile_name)

def get_align_dna_to_alignment_commands(bad_dna_seqs_file, codonalignment_file, 
                                        outdir, aligner):
    geneName = codonalignment_file.split('/')[-1].split('.')[0]
    if aligner == "prank":
        raise Exception("This is a bug! Panaroo supports codon alignment with MAFFT and Clustal only")
    elif aligner == "mafft":
        command = ["mafft",
                   "--add",
                   bad_dna_seqs_file,
                   codonalignment_file,
                   outdir + "aligned_gene_sequences/" + geneName + ".aln.fas"]
    #Note that the MAFFT command must be run with command[:-1] as it writes 
    # to STDOUT by default. Use capture STDOUT when running with subprocess    
    elif aligner == "clustal":
        command = ["clustalo",
                   "--in", bad_dna_seqs_file,
                   "--profile1", codonalignment_file,
                   "--out", outdir + "aligned_gene_sequences/" + geneName + 
                   "aln.fas"
                   ]
    return (command, bad_dna_seqs_file)

def align_sequences(command, outdir, aligner):
    #Avoid running alignments on single-isolate genes
    if command[0] == None:
        return None
    if aligner == "mafft":
        name = str(command[0]).split()[-1].split("/")[-1].split(".")[0]
        stdout, stderr = command[0]()
        with open(outdir + name + ".aln.fas", "w+") as handle:
            handle.write(stdout)
    elif aligner == "clustal":
        try:
            stdout, stderr = command[0]()
        except Bio.Application.ApplicationError as error:
            inputname = str(command[0]).split("-i")[1].split("-t")[0].strip()
            name = inputname.split("/")[-1]
            print(error)
            if "contains 1 sequence, nothing to align" in str(error):
                os.rename(inputname, outdir + name)
            else:
                raise Exception("Clustal failed to run on" + inputname)
    else:
        stdout, stderr = command[0]()
    try:
        os.remove(command[1])
    except FileNotFoundError:
        None
    return True

def realign_dna_sequences(command, outdir, aligner):
    if aligner == "prank":
        raise Exception("This is a bug! Please report it. Panaroo supports " + 
                        "codon alignment with MAFFT and Clustal only")    
    elif aligner == "mafft":
        result = subprocess.Popen(command[0][:-1], stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE)
        mafft_out, mafft_err = result.communicate()
        with open(command[0][-1], 'wb') as outhandle:
            outhandle.write(mafft_out)
    elif aligner == "clustal":
        result = subprocess.Popen(command[0])
    
    #Delete the bad DNA seqs file
    try:
        os.remove(command[1])
    except FileNotFoundError:
        None
    return True
    
def multi_align_sequences(commands, outdir, threads, aligner):
    for command in commands:
        if command == None:
            print(command)
    alignment_results = Parallel(n_jobs=threads, prefer="threads")(
        delayed(align_sequences)(x, outdir, aligner) for x in tqdm(commands)
    )

    return True

def multi_realign_sequences(commands, outdir, threads, aligner):
    
    alignment_results = Parallel(n_jobs=threads, prefer="threads")(
        delayed(realign_dna_sequences)(x, outdir, aligner) for x in tqdm(commands))
    
    return True

def replace_last(string, find, replace):
    reversed = string[::-1]
    replaced = reversed.replace(find[::-1], replace[::-1], 1)
    return replaced[::-1]

def read_sequences(handle):
    with open(handle, 'r') as inhandle:
        sequences = list(SeqIO.parse(inhandle, 'fasta'))
    return sequences

def read_alignment(handle):
    with open(handle, 'r') as inhandle:
        alignment = AlignIO.read(inhandle, 'fasta')
    return alignment

def reverse_translate_sequences(protein_sequence_files, dna_sequence_files, 
                                outdir, temp_directory, aligner, threads):
    #Check that the dna and protein files match up
    for index in range(len(protein_sequence_files)):
        gene_id = protein_sequence_files[index].split('/')[-1].split(".")[0]
        if gene_id == dna_sequence_files[index].split('/')[-1].split(".")[0]:
            continue
        else:
            print(protein_sequence_files[index])
            print(dna_sequence_files[index])
            raise ValueError("DNA and protien sequence IDs do not match!")
    
    #Read in files (multithreaded)
    dna_sequences = Parallel(n_jobs=threads, prefer="threads")(

            delayed(read_sequences)(x) 
            for x in dna_sequence_files)  
    protein_alignments = Parallel(n_jobs=threads, prefer="threads")(
            delayed(read_alignment)(x) 
            for x in protein_sequence_files)
    
    #Check that protein and DNA sequences match, output 
    #Remove DNA sequences that do not match and output to elsewhere, for
    #secondary alignment to sequences aligned at the protein level
    
    clean_dna = []
    clean_proteins = []
    
    reject_dna_files = {}
    
    for index in range(len(dna_sequences)):
        dna = list(dna_sequences[index])
        protein = protein_alignments[index]
        seqids_to_remove = []
        
        #set up sequentially checked QC failure variables for each sequence
        fail_condition_1 = False
        fail_condition_2 = False
        fail_condition_3 = False

        for seq_index in range(len(dna)):
            #Need to take protein without proceeding or trailing gaps
            nogapped_protein_seq = str(protein[seq_index].seq).replace("-", "")
            translated_dna = dna[seq_index].seq.translate()
            
            #fail if the translated sequence isn't the same as the protein
            fail_condition_1 = str(translated_dna).strip("*") != str(nogapped_protein_seq)
            
            #fail if there is a run of > 1 unknown nucleotides
            if fail_condition_1 == False:
                #only test if it hasn't already failed
                fail_condition_2 = "NN" in str(dna[seq_index].seq)
            
            #Fail if the DNA contains degenerate codon, codonalign cannot cope
            if (fail_condition_1 and fail_condition_2) == False:
                #Most expensive test, only test things passing both
                fail_condition_3 = False
                #Such an expensive test, do a cheaper filtering first
                if "N" in str(dna[seq_index].seq):
                    for codon in unambiguous_degenerate_codons.keys():
                        if codon in dna[seq_index].seq:
                            fail_condition_3 = True
            
            if fail_condition_1 or fail_condition_2 or fail_condition_3:
                seqids_to_remove = seqids_to_remove + list(set([dna[seq_index].id, 
                                                                protein[seq_index].id]))
        reject_dna = []        
        #Do the removal if any DNA sequences fail tests
        if (len(seqids_to_remove) > 0):
            clean_nucs = []
            clean_prots = []
            for sequence in dna:
                if sequence.id in seqids_to_remove:
                    reject_dna.append(sequence)
                else:
                    clean_nucs.append(sequence)
            for sequence in protein:
                if sequence.id in seqids_to_remove:
                    continue
                else:
                    clean_prots.append(sequence)
            
            clean_alignment = MultipleSeqAlignment(clean_prots)
            clean_dna.append(clean_nucs)
            clean_proteins.append(clean_alignment)  
            
            gene_name = dna_sequence_files[index].split('/')[-1].split(".")[0]
            reject_outname = temp_directory + gene_name + "_untrans_dna.fasta"
            SeqIO.write(reject_dna, reject_outname, "fasta")
            reject_dna_files[gene_name] = reject_outname
                          
        else:
            clean_dna.append(dna)
            clean_proteins.append(protein)                         
    
    #build codon alignments
    
    #codon_alignments = Parallel(n_jobs=threads, prefer="threads")(
    #        delayed(codonalign.build)
    #        (protein_alignments[index], list(dna_sequences[index])) 
    #        for index in range(len(protein_alignments)))
    

    #do it single threaded because of the need to separate aln missing DNA seq
    
    completed_codon_alignments = {}
    missing_sequences_codon_alignments = {}
    for index in range(len(clean_proteins)):
        gene_name = dna_sequence_files[index].split('/')[-1].split(".")[0]
        try:
            alignment = codonalign.build(clean_proteins[index], clean_dna[index])
            if gene_name in reject_dna_files.keys():
                missing_sequences_codon_alignments[gene_name] = alignment
            else:
                completed_codon_alignments[gene_name] = alignment
        except RuntimeError as e:
            print(e)
            print(index)
            print(protein_sequence_files[index])
            print(dna_sequence_files[index])
        except IndexError as e:
            print(e)
            print(index)
            print(protein_sequence_files[index])
            print(dna_sequence_files[index])
    
    
    #Remove <unknown description> from codon alignments
    for gene in completed_codon_alignments:
        for sequence in completed_codon_alignments[gene]:
            sequence.description = ""
    
    for gene in missing_sequences_codon_alignments:
        for sequence in missing_sequences_codon_alignments[gene]:
            sequence.description = ""
    
    #output successful codon alignments
    write_success_failures = Parallel(n_jobs=threads, prefer="threads")(
            delayed(AlignIO.write)
            (completed_codon_alignments[x], 
             outdir + "aligned_gene_sequences/" + x +".aln.fas", 'fasta')
            for x in completed_codon_alignments)
    
    #output alignments missing some DNA sequences to tmpdir
    
    write_success_failures2 = Parallel(n_jobs=threads, prefer="threads")(
            delayed(AlignIO.write)
            (missing_sequences_codon_alignments[x], 
             temp_directory + x +".aln.fas", 'fasta')
            for x in missing_sequences_codon_alignments)    
    
    print(str(len(missing_sequences_codon_alignments)) + " DNA realignments to perform...")
    
    
    #realign DNA sequences to failed alignments
    
    dna2codons_commands = []
    for gene_name in reject_dna_files:
        command = get_align_dna_to_alignment_commands(reject_dna_files[gene_name], 
                            temp_directory + gene_name + ".aln.fas", 
                            outdir, aligner)
        dna2codons_commands.append(command)
    
    multi_realign_sequences(dna2codons_commands, outdir + "aligned_gene_sequences/",
                              threads, aligner)
            
    
    all_alignments = os.listdir(outdir + "aligned_gene_sequences/")
    
    return all_alignments

def write_alignment_header(alignment_list, outdir, filename):
    out_entries = []
    # Set the tracking variables for gene positions
    gene_start = 1
    gene_end = 0
    for gene in alignment_list:
        # Get length and name from one sequence in the alignment
        # Set variables that need to be set pre-output
        gene_end += gene[2]
        gene_name = gene[0]
        # Create the 3 line feature entry
        gene_entry1 = (
            "FT   feature         " + str(gene_start) + ".." + str(gene_end) + "\n"
        )
        gene_entry2 = "FT                   /label=" + gene_name + "\n"
        gene_entry3 = "FT                   /locus_tag=" + gene_name + "\n"
        gene_entry = gene_entry1 + gene_entry2 + gene_entry3
        # Add it to the output list
        out_entries.append(gene_entry)
        # Alter the post-output variables
        gene_start += gene[2]
    # Create the header and footer
    header = (
        "ID   Genome standard; DNA; PRO; 1234 BP.\nXX\nFH   Key"
        + "             Location/Qualifiers\nFH\n"
    )
    footer = (
        "XX\nSQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T;" + " 0 other;\n//\n"
    )
    # open file and output
    with open(outdir + filename, "w+") as outhandle:
        outhandle.write(header)
        for entry in out_entries:
            outhandle.write(entry)
        outhandle.write(footer)

    return True
