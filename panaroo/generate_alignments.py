import os
import subprocess
import sys
import re

from joblib import Parallel, delayed
from tqdm import tqdm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.Align.Applications import PrankCommandline
from Bio.Align.Applications import MafftCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
import Bio.Application

import numpy as np


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
    # Counter for the number of sequences to
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


def get_alignment_commands(fastafile_name, outdir, aligner, threads):
    geneName = fastafile_name.split("/")[-1].split(".")[0]
    if aligner == "prank":
        command = PrankCommandline(d=fastafile_name, o=geneName, f=8, codon=True)
    elif threads > 3:
        if aligner == "mafft":
            command = MafftCommandline(
                input=fastafile_name, auto=True, nuc=True, adjustdirection=True
            )
        elif aligner == "clustal":
            command = ClustalOmegaCommandline(
                infile=fastafile_name,
                outfile=outdir + "aligned_gene_sequences/" + geneName + ".aln.fas",
                seqtype="DNA",
            )
    elif threads <= 3:
        if aligner == "mafft":
            command = MafftCommandline(
                input=fastafile_name,
                auto=True,
                thread=threads,
                nuc=True,
                adjustdirection=True,
            )
        elif aligner == "clustal":
            command = ClustalOmegaCommandline(
                infile=fastafile_name,
                outfile=outdir + "aligned_gene_sequences/" + geneName + ".aln.fas",
                seqtype="DNA",
                threads=threads,
            )
    return (command, fastafile_name)


def align_sequences(command, outdir, aligner):
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


def multi_align_sequences(commands, outdir, threads, aligner):

    alignment_results = Parallel(n_jobs=threads, prefer="threads")(
        delayed(align_sequences)(x, outdir, aligner) for x in tqdm(commands)
    )

    return True


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


def is_conserved(seqs, threshold):
    # compare everything to the first sequence
    ref = np.fromstring(str(seqs[0][1]), dtype=np.int8)
    length = float(len(ref))
    for seq in seqs[1:]:
        if (
            np.sum(np.fromstring(str(seq[1].lower()), dtype=np.int8) == ref) / length
            < threshold
        ):
            return False
    return True
